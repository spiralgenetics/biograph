#!/bin/sh

# This is a kludge to use the proper system included boost.  Ideally
# we should compile it as part of the workspace instead.

# Autodetect gcc version
set $(gcc --version | head -1)
GCC_VERSION=$3

if echo $GCC_VERSION | egrep '^4'
then
    # GCC version 4.  Use included extra libs.  Corresponding header
    # files are installed in /usr/local.

    cat > BUILD <<EOF

alias(name="boost", actual="@//vendor/extralibs:boost",
visibility=["//visibility:public"],
)
alias(name="python", actual="@//vendor/extralibs:boost_python",
visibility=["//visibility:public"],
)

EOF

elif echo $GCC_VERSION | egrep '^5'
then
    # GCC version 5.  Use system libs.

    ln -s /usr/lib/x86_64-linux-gnu/libboost_*.so .

    cat > BUILD <<EOF

cc_library(name="boost", srcs=[
"libboost_filesystem.so",
"libboost_program_options.so",
"libboost_regex.so",
"libboost_system.so",
"libboost_thread.so",
], linkstatic=1,
visibility=["//visibility:public"],
linkopts=["-lpthread"]
)

cc_library(name="python", srcs=[
"libboost_python.so",
], linkstatic=1,
visibility=["//visibility:public"],
linkopts=["-lpython3.5"]
)


EOF
else
    echo "Unable to autodetect GCC Version: $GCC_VERSION"
exit 1
fi

