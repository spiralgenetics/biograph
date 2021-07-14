cc_library(
    name = "urlencode",
    srcs = ["urlencode.cpp"],
    hdrs = ["urlencode.h"],
    visibility = ["//visibility:public"],
    deps = [
        "//modules/io",
    ],
)

cc_library(
    name = "web",
    srcs = [
        "chunked_encoding.cpp",
        "forward_request_handler.cpp",
        "httpclient.cpp",
        "httpserver.cpp",
        "restful.cpp",
        "url_query.cpp",
        "urlencode.cpp",
    ],
    hdrs = [
        "chunked_encoding.h",
        "couchdb.h",
        "forward_request_handler.h",
        "httpclient.h",
        "httpserver.h",
        "jsontypes.h",
        "rest_exceptions.h",
        "restful.h",
        "url_query.h",
        "urlencode.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "//modules/io",
        "@poco//:net",
        "@poco//:netssl",
    ],
)

cc_library(
    name = "web_ssl",
    srcs = [
        "httpsclient.cpp",
    ],
    hdrs = ["httpsclient.h"],
    visibility = ["//visibility:public"],
    deps = [
        ":web",
        "@poco//:netssl",
    ],
)

cc_test(
    name = "query_variables_test",
    srcs = ["query_variables_test.cpp"],
    deps = [
        "//modules/test:gtest_main",
        "//modules/web",
    ],
)

cc_test(
    name = "urlencode_test",
    srcs = ["urlencode_test.cpp"],
    deps = [
        "//modules/test:gtest_main",
        "//modules/web",
    ],
)
