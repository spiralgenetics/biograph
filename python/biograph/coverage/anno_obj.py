"""
Generic annotation object for organizing
"""

# Template for generating annotation storage classes.
_class_template = """\

class {typename}(object):
    '{typename}()'

    __slots__ = {slotlist!r}

    def __init__(self):
        'Create new instance of {typename} with default values'
{set_defaults}

    @staticmethod
    def get_header():
        'Returns a list of FORMAT lines to be included in the header.'
        return {header_lines!r}

    @staticmethod
    def get_format_tags():
        'Returns a list of tags to be included in the FORMAT field.'
        return {format_tags!r}

    def to_samp_dict(self):
        'Returns a dict of tag to value'
        return {samp_dict}
"""

_set_defaults_template = """
        self.{propname} = {default_value}
"""

def make_annotation_base(typename, annotations):
    """
Generates a "namedtuple"-style class supporting coverage annotations.
"annotations" is a list of tuples of (property, tag, info_desc, default value).
For instance:
 ("upstream_span", "US", "Length of best-spanning read upstream", 0)

"""
    set_defaults = []
    header_lines = []
    format_tags = []
    slotlist = []
    samp_dict = ""
    for propname, tag, info_desc, default_value in annotations:
        set_defaults.append(_set_defaults_template.format(propname=propname,
                                                          default_value=default_value))
        header_lines.append(f"##FORMAT=<ID={tag},Number=R,Type=Integer,Description=\"{info_desc}\">")
        format_tags.append(tag)
        slotlist.append(propname)
        if samp_dict:
            samp_dict += ", "
        samp_dict += f"{tag!r}: self.{propname}"
    classdef = _class_template.format(typename=typename,
                                      set_defaults="".join(set_defaults),
                                      header_lines=header_lines,
                                      format_tags=format_tags,
                                      slotlist=slotlist,
                                      samp_dict="{" + samp_dict + "}")
    namespace = dict(__name__=f'annotation_base_{typename}')
    exec(classdef, namespace) # pylint: disable=exec-used
    result = namespace[typename]
    result._source = classdef
    return result
