#!/usr/bin/env python
"""Generate the parts of Doc/internals_{JP,EN}.typ that are derived from the source.

  api_reference.typ : every public function of libs.plibs / libs.flibs (the __all__ of each
                      submodule) with its arguments and the first sentence of its docstring
                      (@brief if present)
  callgraph.typ     : doxygen-style call graphs drawn with the Typst package diagraph
                      (Graphviz compiled to WASM, so no `dot` binary is needed):
                      the module dependency graph and, for every option of main.py, the
                      calls from main.py down through plibs and flibs to the Fortran kernels

The sources are parsed with ast, so the Fortran library need not be built.  Rerun after
changing the code:

    python Doc/gen_internals.py
"""
import ast
import collections
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parent.parent
PACKAGES = ('plibs', 'flibs')
API_OUT = ROOT / 'Doc' / 'api_reference.typ'
GRAPH_OUT = ROOT / 'Doc' / 'callgraph.typ'
MAXLEN = 220
DIAGRAPH = '#import "@preview/diagraph:0.3.6": raw-render'
GRAPH_DEPTH = 5      # levels drawn below the option node: run_* of main.py -> plibs driver -> ...
GRAPH_MAX_NODES = 40
# argument-marshalling helpers of flibs._loader: called by every wrapper, they would only add noise
LOADER_HELPERS = {'i64', 'dbl', 'c_array', 'f_array', 'cplx', 'c_bool', 'boolean'}
COLORS = {'main': '#e8e8e8', 'plibs': '#dbe9f7', 'flibs': '#dcf0dc', 'fortran': '#fbe3c8'}


def summary(doc: str) -> str:
    m = re.search(r'@brief\s+(.*?)(?:\n\s*\n|\n\s*@|$)', doc, re.S)
    text = ' '.join((m.group(1) if m else doc.split('\n\n')[0]).split())
    m = re.match(r'(.+?[.!?])(?:\s+[A-Z(]|$)', text)
    text = m.group(1) if m else text
    return text if len(text) <= MAXLEN else text[:MAXLEN].rsplit(' ', 1)[0] + ' …'


def arg_names(fn: ast.FunctionDef) -> str:
    a = fn.args
    names = [x.arg for x in a.posonlyargs + a.args]
    if a.vararg:
        names.append('*' + a.vararg.arg)
    names += [x.arg for x in a.kwonlyargs]
    if a.kwarg:
        names.append('**' + a.kwarg.arg)
    return ', '.join(names)


def typst_str(s: str) -> str:
    return '"' + s.replace('\\', '\\\\').replace('"', '\\"') + '"'


# ------------------------------------------------------------------ source model

class Module:
    """One submodule of plibs/flibs: its functions, exports and import aliases."""

    def __init__(self, pkg: str, path: pathlib.Path):
        self.pkg, self.name = pkg, path.stem
        self.key = f'{pkg}.{path.stem}'
        self.tree = ast.parse(path.read_text())
        self.doc = ast.get_docstring(self.tree) or ''
        self.funcs = {n.name: n for n in self.tree.body if isinstance(n, ast.FunctionDef)}
        self.exported = None
        for node in self.tree.body:
            if isinstance(node, ast.Assign) and getattr(node.targets[0], 'id', '') == '__all__':
                self.exported = [e.value for e in node.value.elts]
        self.aliases = {}   # local name -> ('pkg', pkg) | ('mod', key) | ('func', key, name) | ('fortran',)
        for node in ast.walk(self.tree):
            if isinstance(node, ast.Import):
                for a in node.names:
                    self._alias(a.asname or a.name, a.name)
            elif isinstance(node, ast.ImportFrom):
                base = node.module or ''
                if node.level:
                    base = f'libs.{pkg}.{base}' if base else f'libs.{pkg}'
                for a in node.names:
                    full = f'{base}.{a.name}'
                    if not self._alias(a.asname or a.name, full):
                        m = re.fullmatch(r'libs\.(plibs|flibs)\.(_\w+)', base)
                        if m:
                            self.aliases[a.asname or a.name] = ('func', f'{m[1]}.{m[2]}', a.name)

    def _alias(self, local: str, full: str) -> bool:
        m = re.fullmatch(r'libs\.(plibs|flibs)', full)
        if m:
            self.aliases[local] = ('pkg', m[1])
            return True
        m = re.fullmatch(r'libs\.(plibs|flibs)\.(_\w+)', full)
        if m:
            self.aliases[local] = ('mod', f'{m[1]}.{m[2]}')
            return True
        return False

    def public(self) -> list:
        return [f for n, f in self.funcs.items() if not n.startswith('_')
                and (self.exported is None or n in self.exported)]


MAIN_KEY = 'main.main'


def load_modules() -> dict:
    """The plibs/flibs submodules, plus main.py (key MAIN_KEY) for its run_* functions."""
    mods = {MAIN_KEY: Module('main', ROOT / 'main.py')}
    for pkg in PACKAGES:
        for path in sorted((ROOT / 'libs' / pkg).glob('_*.py')):
            if path.name != '__init__.py':
                m = Module(pkg, path)
                mods[m.key] = m
    return mods


def fortran_sources() -> dict:
    """Fortran symbol -> source file, from bind(C, name=...) and subroutine/function names."""
    where = {}
    for f in sorted((ROOT / 'libs' / 'src').glob('*.f90')):
        text = f.read_text(errors='replace')
        for m in re.finditer(r"bind\s*\(\s*C\s*,\s*name\s*=\s*['\"](\w+)['\"]", text, re.I):
            where.setdefault(m[1].lower(), f.name)
        for m in re.finditer(r'^\s*(?:\w+\s+)*(?:subroutine|function)\s+(\w+)', text, re.I | re.M):
            where.setdefault(m[1].lower(), f.name)
    return where


class CallGraph:
    def __init__(self, mods: dict):
        self.mods = mods
        self.pkg_names = {pkg: {} for pkg in PACKAGES}   # exported name -> module key
        for m in mods.values():
            if m.pkg in PACKAGES:
                for f in m.public():
                    self.pkg_names[m.pkg].setdefault(f.name, m.key)
        self.fsrc = fortran_sources()
        self.edges = collections.defaultdict(list)          # node -> ordered callees
        for m in mods.values():
            for name, fn in m.funcs.items():
                node = f'{m.key}.{name}'
                for callee in self._callees(m, fn):
                    if callee != node and callee not in self.edges[node]:
                        self.edges[node].append(callee)

    def resolve(self, m: Module, func: ast.expr) -> str | None:
        if isinstance(func, ast.Name):
            if func.id in m.funcs:
                return f'{m.key}.{func.id}'
            a = m.aliases.get(func.id)
            if a and a[0] == 'func' and a[1] in self.mods and a[2] in self.mods[a[1]].funcs:
                return f'{a[1]}.{a[2]}'
            return None
        if not isinstance(func, ast.Attribute):
            return None
        val = func.value
        # _loader._lib.NAME(...) is a call into the Fortran library
        if isinstance(val, ast.Attribute) and val.attr == '_lib':
            return 'fortran.' + func.attr.lower()
        if isinstance(val, ast.Name):
            if m.key == 'flibs._loader' and val.id == '_lib':
                return 'fortran.' + func.attr.lower()
            a = m.aliases.get(val.id)
            if a is None:
                return None
            if a[0] == 'pkg':
                key = self.pkg_names[a[1]].get(func.attr)
                return f'{key}.{func.attr}' if key else None
            if a[0] == 'mod' and a[1] in self.mods and func.attr in self.mods[a[1]].funcs:
                if a[1] == 'flibs._loader' and func.attr in LOADER_HELPERS:
                    return None
                return f'{a[1]}.{func.attr}'
        return None

    def _callees(self, m: Module, fn: ast.FunctionDef):
        for node in ast.walk(fn):
            if isinstance(node, ast.Call):
                target = self.resolve(m, node.func)
                if target:
                    yield target

    def kind(self, node: str) -> str:
        return node.split('.', 1)[0]

    def label(self, node: str) -> str:
        if node.startswith('fortran.'):
            sym = node.split('.', 1)[1]
            src = self.fsrc.get(sym) or self.fsrc.get(sym.rstrip('_'))
            return f'{sym}\\n[{src}]' if src else sym
        if node.startswith(MAIN_KEY + '.'):
            return node.rsplit('.', 1)[1] + '\\nmain.py'
        if node.startswith('main.'):
            return node.split('.', 1)[1]
        pkg, mod, name = node.split('.')
        return f'{name}\\n{pkg}.{mod}'


# ------------------------------------------------------------------ main.py dispatch

def main_dispatch() -> list:
    """[(option value, CalcMode name, [run_* function node])] from the RUNNERS table of main.py."""
    tree = ast.parse((ROOT / 'main.py').read_text())
    values, runners = {}, {}
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == 'CalcMode':
            for st in node.body:
                if isinstance(st, ast.Assign) and isinstance(st.value, ast.Tuple):
                    values[st.targets[0].id] = st.value.elts[0].value
        elif isinstance(node, ast.Assign) and getattr(node.targets[0], 'id', '') == 'RUNNERS':
            for k, v in zip(node.value.keys, node.value.values):
                runners[k.attr] = [f'{MAIN_KEY}.{v.id}']
    return sorted(((values[k], k, v) for k, v in runners.items()), key=lambda t: t[0])


# ------------------------------------------------------------------ dot output

def node_attr(graph: CallGraph, node: str) -> str:
    kind = graph.kind(node)
    shape = 'box' if kind != 'fortran' else 'box, style="rounded,filled"'
    style = '' if kind == 'fortran' else ', style=filled'
    return f'"{node}" [label="{graph.label(node)}", shape={shape}{style}, fillcolor="{COLORS[kind]}"];'


def option_dot(graph: CallGraph, value: int, name: str, drivers: list) -> str:
    root = f'main.option={value} {name}'
    order, depth, edges = [root], {root: 0}, []
    queue = collections.deque()
    for d in drivers:
        edges.append((root, d))
        if d not in depth:
            depth[d] = 1
            order.append(d)
            queue.append(d)
    truncated = set()
    while queue:
        node = queue.popleft()
        for callee in graph.edges.get(node, []):
            if callee not in depth:
                if depth[node] >= GRAPH_DEPTH or len(order) >= GRAPH_MAX_NODES:
                    truncated.add(node)
                    continue
                depth[callee] = depth[node] + 1
                order.append(callee)
                queue.append(callee)
            edges.append((node, callee))
    lines = ['digraph {', '  rankdir=LR; nodesep=0.12; ranksep=0.35;',
             '  node [fontname="Helvetica", fontsize=9, margin="0.06,0.03", height=0.3];',
             '  edge [arrowsize=0.5, color="#555555"];']
    lines += ['  ' + node_attr(graph, n) for n in order]
    for n in sorted(truncated):
        lines.append(f'  "{n}__more" [label="…", shape=plaintext];')
        lines.append(f'  "{n}" -> "{n}__more" [style=dashed];')
    lines += [f'  "{a}" -> "{b}";' for a, b in edges]
    lines.append('}')
    return '\n'.join(lines)


def module_dot(graph: CallGraph, layer: str) -> str:
    """Module-level dependencies, edge width ~ number of distinct function-to-function calls.
    layer='plibs': the plibs submodules among themselves;
    layer='kernels': plibs -> flibs -> Fortran source files (calls inside plibs left out)."""
    def mod_of(node: str) -> str:
        if node.startswith('fortran.'):
            sym = node.split('.', 1)[1]
            return 'fortran.' + (graph.fsrc.get(sym) or graph.fsrc.get(sym.rstrip('_')) or '?')
        return node.rsplit('.', 1)[0]
    weight = collections.Counter()
    for a, callees in graph.edges.items():
        for b in callees:
            ma, mb = mod_of(a), mod_of(b)
            if ma == MAIN_KEY:
                continue
            both_plibs = ma.startswith('plibs.') and mb.startswith('plibs.')
            if ma != mb and both_plibs == (layer == 'plibs'):
                weight[(ma, mb)] += 1
    mods = sorted({m for e in weight for m in e})
    lines = ['digraph {', f'  rankdir={"TB" if layer == "plibs" else "LR"}; nodesep=0.12; ranksep=0.5;',
             '  node [fontname="Helvetica", fontsize=9, margin="0.06,0.03", height=0.28, style=filled];',
             '  edge [arrowsize=0.45, color="#777777"];']
    for kind in ('plibs', 'flibs', 'fortran'):
        members = [m for m in mods if m.startswith(kind + '.')]
        if not members:
            continue
        title = 'libs/src (Fortran)' if kind == 'fortran' else kind
        lines.append(f'  subgraph cluster_{kind} {{ label="{title}"; fontname="Helvetica"; fontsize=10; color="#bbbbbb";')
        for m in members:
            shape = 'box, style="rounded,filled"' if kind == 'fortran' else 'box'
            lines.append(f'    "{m}" [label="{m.split(".", 1)[1]}", shape={shape}, fillcolor="{COLORS[kind]}"];')
        lines.append('  }')
    for (a, b), w in sorted(weight.items()):
        lines.append(f'  "{a}" -> "{b}" [penwidth={min(0.5 + 0.25 * w, 3):.2f}];')
    lines.append('}')
    return '\n'.join(lines)


def render(dot: str) -> str:
    fence = '`' * 6
    return f'raw-render({fence}\n{dot}\n{fence})'


# Typst helpers written at the top of callgraph.typ.  fit() draws a graph at its natural size
# and scales it down (never up) to the column width and to GRAPH_MAX_HEIGHT.  fit-here() does the
# same for a heading + graph kept together, and if that would not fit in the space left on the
# current page after the previous graph it shrinks the graph into that space as long as it keeps >= GRAPH_MIN_SHRINK of the
# fitted size (the labels stay readable); otherwise the pair moves to the next page.  This keeps pages
# from ending half empty.  The page margins are those of internals_{JP,EN}.typ (2 cm).
FIT = """#let fit(g, max-height: GRAPH_MAX_HEIGHT) = layout(size => {
  let s = measure(g)
  let k = calc.min(1.0, size.width / s.width, max-height / s.height)
  align(center, scale(k * 100%, reflow: true, g))
})

#let fit-here(head, g, max-height: GRAPH_MAX_HEIGHT, min-shrink: GRAPH_MIN_SHRINK) = context {
  // the space left is measured from the end of the previous graph: a probe placed before this
  // block would travel with it to the next page, the end of the previous block stays put
  let prev = query(metadata.where(value: "fit-end").before(here()))
  let s = measure(g)
  let k = calc.min(1.0, (page.width - 4cm) / s.width, max-height / s.height)
  if prev.len() > 0 {
    let left = page.height - 2cm - prev.last().location().position().y - 1.6cm
    if s.height * k > left and s.height * k * min-shrink <= left {
      k = left / s.height
    }
  }
  block(breakable: false, width: 100%)[
    #head
    #align(center, scale(k * 100%, reflow: true, g))
    #metadata("fit-end")
  ]
}"""
GRAPH_MAX_HEIGHT = '18cm'   # about 3/4 of the 25.7 cm text height of A4 with 2 cm margins
GRAPH_MIN_SHRINK = '0.65'  # labels may become at most 30% smaller than in the page-fitted graph


def write_callgraph(graph: CallGraph) -> None:
    """callgraph.typ defines module-graph-plibs / module-graph-kernels (one figure each) and
    option-graphs (one level-3 heading + figure per option, kept together by fit-here);
    the including document places them under its own headings."""
    lines = ['// Generated by Doc/gen_internals.py; do not edit by hand.', DIAGRAPH, '',
             FIT.replace('GRAPH_MAX_HEIGHT', GRAPH_MAX_HEIGHT).replace('GRAPH_MIN_SHRINK', GRAPH_MIN_SHRINK), '',
             f'#let module-graph-plibs = fit({render(module_dot(graph, "plibs"))})', '',
             f'#let module-graph-kernels = fit({render(module_dot(graph, "kernels"))})', '',
             '#let option-graphs = [']
    for value, name, drivers in main_dispatch():
        lines += [f'#fit-here([=== `option={value}` (`{name}`)], {render(option_dot(graph, value, name, drivers))})', '']
    lines.append(']')
    GRAPH_OUT.write_text('\n'.join(lines))
    print(f'wrote {GRAPH_OUT}')


def write_api_reference(mods: dict) -> None:
    lines = ['// Generated by Doc/gen_internals.py from the docstrings; do not edit by hand.',
             '#let fn(sig, doc) = block(breakable: false, below: 0.9em)[',
             '  #text(size: 8.5pt, weight: "bold", raw(sig)) \\',
             '  #pad(left: 1.2em, text(size: 9.5pt, doc))',
             ']', '']
    for pkg in PACKAGES:
        lines.append(f'== `libs.{pkg}`')
        lines.append('')
        for m in (m for m in mods.values() if m.pkg == pkg):
            fns = m.public()
            if not fns:
                continue
            lines.append(f'=== `{m.key}`')
            lines.append('')
            if m.doc:
                lines.append(f'#text(size: 9.5pt)[#{typst_str(summary(m.doc))}]')
                lines.append('')
            for fn in fns:
                sig = typst_str(f'{fn.name}({arg_names(fn)})')
                doc = typst_str(summary(ast.get_docstring(fn) or ''))
                lines.append(f'#fn({sig}, {doc})')
            lines.append('')
    API_OUT.write_text('\n'.join(lines))
    print(f'wrote {API_OUT}')


def main() -> None:
    mods = load_modules()
    write_api_reference(mods)
    write_callgraph(CallGraph(mods))


if __name__ == '__main__':
    main()
