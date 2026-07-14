import sys, io, nbformat, contextlib, traceback

def execute(path):
    nb = nbformat.read(path, as_version=4)
    ns = {"__name__": "__main__"}
    count = 0
    for cell in nb.cells:
        if cell.cell_type != "code":
            continue
        count += 1
        cell.execution_count = count
        cell.outputs = []
        buf = io.StringIO()
        err = None
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            try:
                exec(compile(cell.source, f"<cell {count}>", "exec"), ns)
            except Exception:
                err = traceback.format_exc()
        text = buf.getvalue()
        if text:
            cell.outputs.append(nbformat.v4.new_output("stream", name="stdout", text=text))
        if err:
            cell.outputs.append(nbformat.v4.new_output(
                "error", ename="Error", evalue="cell failed",
                traceback=err.splitlines()))
            nbformat.write(nb, path)
            return err
    nbformat.write(nb, path)
    return None

import sys
e = execute(sys.argv[1] if len(sys.argv)>1 else "notebook.ipynb")
if e:
    print("CELL_ERROR:\n", e[-1500:]); sys.exit(1)
print("EXECUTED_OK")
