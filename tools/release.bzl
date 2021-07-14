def _spiral_release_impl(ctx):
    versioned_base = ctx.attr.package_name + "-" + ctx.attr.version
    build_script = """
set -e
TEMP_DIR=$(mktemp -d)
function cleanup {
    rm -r ${TEMP_DIR}
}
trap cleanup EXIT

TAR_DIR='__VERSIONED_BASE__'
TAR_DIR_PATH=${TEMP_DIR}/${TAR_DIR}

mkdir ${TAR_DIR_PATH}
"""
    build_script = build_script.replace("__VERSIONED_BASE__", versioned_base)

    if ctx.attr.file_dirs:
        if len(ctx.attr.file_dirs) != len(ctx.attr.files):
            fail("Length of file_dirs must match length of files", "file_dirs")
        for file_dir in ctx.attr.file_dirs:
            mkdir_cmd = "mkdir -p ${TAR_DIR_PATH}/'" + file_dir + "'\n"
            build_script = build_script + mkdir_cmd

    for input_index in range(len(ctx.attr.files)):
        input_label = ctx.attr.files[input_index]
        for file_info in input_label.files.to_list():
            install_dir = "."
            if ctx.attr.file_dirs:
                install_dir = ctx.attr.file_dirs[input_index]
            copy_cmd = "cp '" + file_info.path + "' ${TAR_DIR_PATH}/'" + install_dir + "/'\n"
            build_script = build_script + copy_cmd
    for link_src, link_dest in ctx.attr.links.items():
        link_cmd = "ln -s '" + link_src + "' ${TAR_DIR_PATH}/'" + link_dest + "'\n"
        build_script = build_script + link_cmd
    build_script = build_script + """
# fix up permissions and strip binaries
find ${TEMP_DIR} -type f -perm /a+x -exec chmod 0755 {} \; -exec strip {} \; > /dev/null 2>&1
find ${TEMP_DIR} -type f ! -perm /a+x -exec chmod 0644 {} \;
find ${TEMP_DIR} -type d -exec chmod 0755 {} \;
"""
    if ctx.attr.pre_tar:
        build_script = build_script + "pushd ${TEMP_DIR} > /dev/null\n" + ctx.attr.pre_tar + "\npopd > /dev/null\n"

    build_script = build_script + "\ntar -zcf __OUT_TARFILE__ -C ${TEMP_DIR} ${TAR_DIR}\n"

    build_script = build_script.replace("__OUT_TARFILE__", ctx.outputs.out_tarfile.path)
    input_list = []
    for target in ctx.attr.files:
        for file_info in target.files.to_list():
            input_list.append(file_info)
    ctx.actions.run_shell(
        inputs=input_list,
        outputs=[ctx.outputs.out_tarfile],
        progress_message = "Packaging " + ctx.outputs.out_tarfile.short_path,
        command = build_script)
    return struct(runfiles=ctx.runfiles(files=[ctx.outputs.out_tarfile]))

_spiral_release_rule = rule(
    attrs = {
        "package_name": attr.string(),
        "version": attr.string(),
        "files": attr.label_list(allow_files = True),
        "file_dirs": attr.string_list(),
        "links": attr.string_dict(),
        "out_tarfile": attr.output(),
        "pre_tar": attr.string()
    },
    implementation = _spiral_release_impl,
)

def spiral_release(name, version, package_name = None, pre_tar = None, files = [], file_dirs = [],
                   links = {}, visibility=[], tags=[]):
    """Prepares a release tarball.

version:  Version to tag the tarball with.

package_name:  Base name of the package to use in naming the tarball and
directory inside.

files: A list of targets to include in the package.

file_dirs: If present, this should be a list of relative directories.  Each directory
in this list should correspond to a target in "files".  The directory specifies
the subdirectory within the tarball that files from this target should be put in.

links: A dictionary of symbolic links that should be generated in a tarball, in
the format {"real_file": "alias"}.

"""
    if not package_name:
        package_name = name
    tarfile_name = package_name + "-" + version + ".tgz"
    _spiral_release_rule(name = name,
                         package_name = package_name,
                         version = version,
                         files = files,
                         file_dirs = file_dirs,
                         links = links,
                         out_tarfile = tarfile_name,
                         visibility = visibility,
                         tags = tags,
                         pre_tar = pre_tar)
