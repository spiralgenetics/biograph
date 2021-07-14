# Bazel repository that autodetects whether we should use
# our custom boost libraries (for gcc 4.8) or stock system
# boost libraries (gcc 5).

def _impl(repository_ctx):
    repository_ctx.execute([
	repository_ctx.path(Label("//tools:boost_repository.sh"))]
    )

# "version" attribute used to force regeneration of the repository
# if something changes.
spiral_boost_repository = repository_rule(
    implementation=_impl,
	attrs={"version": attr.int()})
    
