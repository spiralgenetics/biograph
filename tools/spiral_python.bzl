load("@subpar//:subpar.bzl", "par_binary")
load(
    "@rules_python//python:python.bzl",
      "py_binary", "py_library", "py_test",
)

load("@tool_requirements//:requirements.bzl",
     _tool_requirement="requirement")
load("python_repository.bzl", "PYTHON_VERSIONS")

def spiral_py_library(name, srcs=[], data=[], deps=[], lint_deps=[], **kwargs):
    py_library(name=name, data=data, srcs=srcs,
               deps=deps + ["//python:python_import_base"], **kwargs)
    _lint_srcs(name=name, data=data, srcs=srcs, deps=deps + lint_deps)

# A standalone python executable that doesn't need to be packaged up.

def spiral_py_binary(name, srcs, data=[], deps=[], **kwargs):
    py_binary(name=name, data=data, srcs=srcs, python_version="PY3",
              deps=deps + ["//python:python_import_base"], **kwargs)
    _lint_srcs(name=name, data=data, srcs=srcs, deps=deps)

# A python executable packaged into a python archive (PAR) file.


def spiral_par_binary(name, srcs, data=[], deps=[], **kwargs):
    par_binary(name=name, data=data, srcs=srcs,
               deps=deps + ["//python:python_import_base"], **kwargs)
    _lint_srcs(name=name, data=data, srcs=srcs, deps=deps)


def spiral_py_test(name, srcs, data=[], deps=[], tags=[], **kwargs):
    py_test(name=name, data=data, srcs=srcs,
            deps=deps + ["//python:python_import_base"],
            tags=tags + ["noasan", "py_test"], **kwargs)
    _lint_srcs(name=name, data=data, srcs=srcs, deps=deps)


def _lint_srcs(name, data, srcs, deps):
    if not srcs:
        return
    pkgs = []
    found_init = None
    for src in srcs:
        if ":" not in src:
            src = "//" + native.package_name() + ":" + src
        src = Label(src)
        pkg_name = src.package.replace("/", ".")
        if pkg_name[:7] == "python.":
            pkg_name = pkg_name[7:]
        if src.name == "__init__.py":
            pkgs.append(pkg_name)
            found_init = pkg_name
        elif src.name[-3:] != ".py":
            fail("Expected source " + src.name + " to end in .py")
        else:
            pkgs.append(pkg_name + "." + src.name[:-3])

    if found_init:
        _lint_packages(name=found_init, pkgs=pkgs, data=data, deps=deps)
        return

    for src in srcs:
        _lint_file(src=src, data=data, deps=deps)


def _lint_file(src, data, deps):
    name = src.replace(".", "_").replace("/", "_").replace(":", "_")
    if name[0:2] == "__":
        name = name[2:]
    _lint_internal(name, src=src, data=data, deps=deps)


def _lint_packages(name, pkgs, data, deps):
    name = name.replace(".", "_").replace("/", "_")
    _lint_internal(name, pkgs=pkgs, data=data, deps=deps)


def _lint_internal(name, src=None, pkgs=None, data=None, deps=None):
    cmd = """
if $(location """ + name + """_lint_runner) --rcfile=$(location //tools:pylintrc)"""
    src_deps = []
    if src:
        if pkgs:
            fail("_lint_internal must not get both src and pkg")
        cmd = cmd + " $(location " + src + ")"
        src_deps = [src]
    elif pkgs:
        for pkg in pkgs:
            cmd = cmd + " " + pkg
    else:
        fail("_lint_internal must get either src or pkg")

    cmd = cmd + """ > $(location """ + name + """_lint_output.txt) 2>&1
then
  exit_status=0
else
  exit_status=$$?
fi
if grep useless-suppression $(location """ + name + """_lint_output.txt)
then
  exit_status=1
fi
if grep 'Problem importing module' $(location """ + name + """_lint_output.txt)
then
  exit_status=1
fi
echo $${exit_status} > $(location """ + name + """_lint.status) 2>&1
cat > $(location """ + name + """_lint_test.sh) << 'EOF'
cat $$1
exit `cat $$2`
EOF

exit 0
"""

    py_binary(name=name + "_lint_runner",
              srcs=["//tools:run_pylint.py"],
              main="//tools:run_pylint.py",
              python_version="PY3",
              deps=depset(deps + [_tool_requirement("pylint"), "//python:python_import_base"]),
              )
    # print(depset(deps + TOOL_REQUIREMENTS +
    #           ["//python:python_import_base"]))
    # Yuck; sh_test isn't flexible enough to do what we want, and
    # genrule doesn't let us mark it as a test.
    native.genrule(name=name + "_lint_output",
                   outs=[name + "_lint_output.txt", name + "_lint.status",
                         name + "_lint_test.sh"],
                   tools=[name + "_lint_runner"],
                   srcs=src_deps + ["//tools:pylintrc"], cmd=cmd)
    native.sh_test(name=name + "_lint_test",
                   srcs=[name + "_lint_test.sh"],
                   args=["$(location " + name + "_lint_output.txt)",
                         "$(location " + name + "_lint.status)"],
                   data=[name + "_lint_output.txt", name + "_lint.status"], tags=["pylint"])


def spiral_py_lint(name, srcs, data=[], deps=[]):
    _lint_srcs(name, data, srcs, deps)


def spiral_py_sdist(name, out, deps, packages, scripts=[], bins=[], requirements=None, readme=None, verbose=False, visibility=[]):
    cmd = "$(location :" + name + "_sdist_maker)"
    genrule_srcs = list(deps)
    for script in scripts:
        cmd = cmd + " --script $(location " + script + ")"
        genrule_srcs.append(script)
    if readme:
        cmd = cmd + " --readme $(location " + readme + ")"
        genrule_srcs.append(readme)
    for bin in bins:
        cmd = cmd + " --bin $(location " + bin + ")"
        genrule_srcs.append(bin)
    cmd = cmd + " --out $(location " + out + ")"
    if requirements:
        cmd = cmd + " --requirements $(location " + requirements + ")"
        genrule_srcs.append(requirements)
    if verbose:
        cmd = cmd + " --verbose"
    for pkg in packages:
        cmd = cmd + " " + pkg

    # We can't use the same make_sdist target for all packages, since
    # make_sdist needs to depend on the packages
    spiral_py_binary(name=name + "_sdist_maker",
                     srcs=["//tools:make_sdist.py"],
                     main="//tools:make_sdist.py",
                     deps=deps + [_tool_requirement("setuptools")])
    native.genrule(name=name, outs=[out], srcs=genrule_srcs, cmd=cmd, tools=[
                   ":" + name + "_sdist_maker"], visibility=visibility)

PythonVersion=provider(
    doc="Python version specifciation",
    fields={
        "version": "python version number"
    }
)

def _spiral_python_version_spec_impl(ctx):
    return PythonVersion(version=ctx.build_setting_value)

spiral_python_version_spec = rule(
    implementation=_spiral_python_version_spec_impl,
    build_setting = config.string()
)

def _multiversion_transition_impl(settings, attr):
    return {version: {"//tools:python_version": version} for version in PYTHON_VERSIONS}

_multiversion_transition = transition(
    implementation = _multiversion_transition_impl,
    inputs = [],
    outputs = ["//tools:python_version"])

def _multiversion_module_impl(ctx):
    versioned_deps = ctx.split_attr.module
    output_list = []
    for version, dep in versioned_deps.items():
        out_name = ctx.actions.declare_file(ctx.attr.name + "_" + version.replace(".","") + ".so")
        output_list.append(out_name)
        for input_file in dep[OutputGroupInfo].interface_library.to_list():
            ctx.actions.symlink(output=out_name, target_file=input_file)
    return [DefaultInfo(runfiles=ctx.runfiles(files=output_list))]

spiral_py_multiversion_module = rule(
    implementation = _multiversion_module_impl,
    attrs = {
        "module": attr.label(cfg=_multiversion_transition),
        "_whitelist_function_transition": attr.label(
                 default = "@bazel_tools//tools/whitelists/function_transition_whitelist"
             ),
    })

def _spiral_python_version_alias_impl(ctx):
    n = 0
    actual = None
    version_setting = ctx.attr.python_version[PythonVersion].version
    if version_setting == 'host':
        actual = ctx.attr.host_actual
    else:
        for k, v in ctx.attr.actuals.items():
            if v == version_setting:
                actual = k
    if not actual:
        fail("Could not find python version " + version_setting)
    if CcInfo in actual:
        return struct(CcInfo=actual[CcInfo])
    if PyRuntimeInfo in actual:
        return struct(PyRuntimeInfo=actual[PyRuntimeInfo])
    fail("Unable to find applicable provider in " + str(actual))

_spiral_python_version_alias = rule(
    implementation=_spiral_python_version_alias_impl,
    attrs = {
        "actuals": attr.label_keyed_string_dict(),
        "host_actual": attr.label(),
        "python_version": attr.label(),
        })

def spiral_python_version_alias(name, actual, host_actual, visibility = []):
    _spiral_python_version_alias(name=name,
                                 actuals={"@spiral_python" + version + actual: version
                                          for version in PYTHON_VERSIONS},
                                 host_actual=host_actual,
                                 python_version="//tools:python_version",
                                 visibility=visibility)

