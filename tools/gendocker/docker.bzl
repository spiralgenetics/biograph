def docker_pip_install_test(name, docker_image, lambda_toy_data, pip_package, model,
                            tags = []):
    native.sh_test(
        name = name,
        srcs = ["install_test.sh"],
        data = [pip_package, "run_install_test.sh"],
        args = [
            name,
            docker_image,
            lambda_toy_data,
            "$(location " + pip_package + ")",
            model
        ],
        tags = tags,
        # Run locally, instead of inside whatever build container we were using.
        local = True,
    )

def docker_pip_install_tests(images, repo, **kwargs):
    for image in images:
        docker_pip_install_test(
            name = image + "_install_test",
            docker_image = repo + ":" + image,
            **kwargs
        )

def centos7_tmpl(name, tag):
    docker_template(name=name,
                    out="Dockerfile-" + name,
                    template="Dockerfile-centos7.tmpl",
                    substitutions={"{TAG}": tag})

def ubuntu_deadsnake_tmpl(name, py_vers, tag):
    docker_template(name=name,
                    out="Dockerfile-" + name,
                    template="Dockerfile-ubuntu-deadsnake.tmpl",
                    substitutions={"{PY_VERS}": py_vers,
                                   "{TAG}": tag})

def ubuntu_tmpl(name, tag):
    docker_template(name=name,
                    out="Dockerfile-" + name,
                    template="Dockerfile-ubuntu.tmpl",
                    substitutions={"{TAG}": tag})

def _docker_template_impl(ctx):
    ctx.actions.expand_template(
        output = ctx.outputs.out,
        template = ctx.file.template,
        substitutions = ctx.attr.substitutions,
    )
    return [DefaultInfo(files = depset([ctx.outputs.out]))]

docker_template = rule(
    attrs = {
        "out": attr.output(mandatory=True),
        "template": attr.label(mandatory=True, allow_single_file=True),
        "substitutions": attr.string_dict(mandatory=True)
    },
    implementation=_docker_template_impl)


def docker_image_builder(*, name, srcs, data, build_images, test_images):
    args = ["--docker-repo-base", "$(DOCKER_REPO_BASE)"]
    for image in build_images:
        args.extend(["--build-image", image])
    for image in test_images:
        args.extend(["--test-image", image])
    native.sh_binary(
        name=name,
        srcs=srcs,
        data=data,
        args=args)

