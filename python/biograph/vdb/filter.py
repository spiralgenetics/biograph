#!/usr/bin/env python3
'''
VDB query filter parser

Syntax is based on common bcftools expressions, adapted for VDB:

  https://samtools.github.io/bcftools/bcftools.html#expressions

The parsing code is based heavily on the fourFn.py example from pyparsing.

TODOs:

bcftools & and | are more complex than AND and OR, since they must simultaneously match
within a single sample. I believe this would require a HAVING clause.

@file notation (consider a vdb query escape character for similar functionality)
binom(), phred()
N_ALT, N_SAMPLES, AC, MAC, AF, MAF, AN, ILEN
N_PASS, F_PASS

'''

from pyparsing import (
    alphas,
    alphanums,
    oneOf,
    pyparsing_common,
    replaceWith,
    CaselessKeyword,
    Combine,
    Forward,
    Keyword,
    Literal,
    ParseException,
    QuotedString,
    Word,
    ZeroOrMore
)

__all__ = ["parser", "ParseException"]

def generate_vdb_parser(): # pylint: disable=too-many-locals, too-many-statements
    ''' The syntax parser definition '''
    # Many more are available in Presto: https://prestodb.io/docs/0.217/functions.html
    Functions = (
        CaselessKeyword('MAX') |
        CaselessKeyword('MIN') |
        CaselessKeyword('AVG') |
        CaselessKeyword('MEAN').setParseAction(replaceWith('GEOMETRIC_MEAN')) |
        # CaselessKeyword('MEDIAN') | # approx_percentile(val, 0.5)
        CaselessKeyword('STDEV').setParseAction(replaceWith('STDDEV')) |
        CaselessKeyword('SUM') |
        CaselessKeyword('STRLEN').setParseAction(replaceWith('LENGTH')) |
        CaselessKeyword('ABS') |
        CaselessKeyword('COUNT')
    )
    Functions.setName('Function')

    # SampleFunctions = (
    #     CaselessKeyword('SMPL_MAX') |
    #     CaselessKeyword('SMPL_MIN') |
    #     CaselessKeyword('SMPL_AVG') |
    #     CaselessKeyword('SMPL_MEAN') |
    #     CaselessKeyword('SMPL_MEDIAN') |
    #     CaselessKeyword('SMPL_STDEV') |
    #     CaselessKeyword('SMPL_SUM') |
    #     CaselessKeyword('sMAX') |
    #     CaselessKeyword('sMIN') |
    #     CaselessKeyword('sAVG') |
    #     CaselessKeyword('sMEAN') |
    #     CaselessKeyword('sMEDIAN') |
    #     CaselessKeyword('sSTDEV') |
    #     CaselessKeyword('sSUM')
    # )

    Number = pyparsing_common.number
    Number.setName('Number')

    # strings must be 'quoted'
    String = (QuotedString('"') | QuotedString("'")).setParseAction(lambda s, l, t: [f"'{t[0]}'"])
    String.setName('String')

    Comparison = Literal("==").setParseAction(replaceWith('=')) | oneOf("= > >= <= < !=")
    Comparison.setName('Comparison')

    # TODO: Keyword('&') | Keyword('|')
    Logical = (
        CaselessKeyword('AND') | Keyword('&&').setParseAction(replaceWith('AND')) |
        CaselessKeyword('OR') | Keyword('||').setParseAction(replaceWith('OR')) |
        Keyword(',').setParseAction(replaceWith('OR')) |
        Keyword('!').setParseAction(replaceWith('NOT'))
    )
    Logical.setName('Logical operator')

    Filter = (CaselessKeyword('filt') | CaselessKeyword('filter')).setParseAction(replaceWith('filt'))
    VCFColumn = (
        CaselessKeyword('chrom') |
        CaselessKeyword('pos') |
        CaselessKeyword('id').setParseAction(replaceWith('varid')) |
        CaselessKeyword('ref') |
        CaselessKeyword('alt') |
        CaselessKeyword('qual') |
        Filter
    )
    VCFColumn.setName('VCF column')

    # All BioGraph FORMAT fields, plus a few common others, CAST() to the proper type.
    FormatInt = oneOf("DP DV GQ LAALTSEQLEN LALANCH LARANCH LAREFSPAN LASCORE NUMASM OV PDP PI RC")
    FormatInt.setParseAction(lambda s, l, t: [f"""CAST(sample['{t[0]}'] AS BIGINT)"""])

    FormatFloat = oneOf("LAALTGC LAREFGC")
    FormatFloat.setParseAction(lambda s, l, t: [f"""CAST(sample['{t[0]}'] AS DOUBLE)"""])

    FormatStr = oneOf("AC AD DC DCC DDC DMO DS DXO EC GT MC MO MP NR PAD PG PL UC UCC UDC UMO US UXO XC XO")
    FormatStr.setParseAction(lambda s, l, t: [f"""sample['{t[0]}']"""])

    FormatField = (FormatInt | FormatFloat | FormatStr)
    FormatField.setName('FORMAT field')

    # All INFO fields, CAST() to the proper type.
    InfoInt = oneOf("SVLEN END")
    InfoInt.setParseAction(lambda s, l, t: [f"""CAST(info['{t[0]}'] AS BIGINT)"""])

    InfoStr = oneOf("SVTYPE")
    InfoStr.setParseAction(lambda s, l, t: [f"""info['{t[0]}']"""])

    InfoField = (InfoInt | InfoStr)
    InfoField.setName('INFO field')

    VCFLookup = (
        Combine(CaselessKeyword('info') + Literal('/') + Word(alphas)).setParseAction(lambda s, l, t: [f"""info['{t[0].split("/")[1]}']"""]) |
        Combine((CaselessKeyword('fmt') | CaselessKeyword('format')) + Literal('/') + FormatField).setParseAction(lambda s, l, t: [f"""{t[0].split("/")[1]}"""])
    )
    VCFLookup.setName('VCF lookup')

    Genotype = Combine(oneOf("0 1 2 .") + oneOf("/ |") + oneOf("0 1 2 .")).setParseAction(lambda s, l, t: [f"'{t[0]}'"])

    # In case you forgot the quotes around your term
    QuoteFix = ((CaselessKeyword('chrom') | Filter) + Comparison + (Word(alphanums))).setParseAction(lambda s, l, t: [f"{t[0]} {t[1]} '{t[2]}'"])
    # Filters are specified in VCF (1-based) coordinates, but stored 0-based.
    PosFix = (CaselessKeyword('pos') + Comparison + Number).setParseAction(lambda s, l, t: [f"{t[0]} {t[1]} {t[2] - 1}"])

    FixUp = (QuoteFix | PosFix)

    # Other VDB columns
    OtherColumns = (
        CaselessKeyword("spans") |
        CaselessKeyword("reflen") |
        CaselessKeyword("varend") |
        CaselessKeyword("varid") |
        CaselessKeyword("checkpoint") |
        CaselessKeyword("study_name") |
        CaselessKeyword("sample_name") |
        CaselessKeyword("aid")
    )

    Null = (Literal('"."') | Literal("'.'")).setParseAction(replaceWith('NULL'))
    IsNull = Literal("=").setParseAction(replaceWith('IS')) + Null

    Field = FormatField | InfoField | VCFColumn | VCFLookup | OtherColumns
    Field.setName('Field')

    ArithOp = oneOf("+ * - /")
    ArithOp.setName('Arithmetic operator')

    Lpar, Rpar = map(Literal, "()")

    # Forward() is a placeholder for a recursive token pattern
    Expression = Forward()
    Expression.setName('Expression')

    ExpressionList = Expression + ZeroOrMore(Literal(',') + Expression).setName('Expression list')
    FunctionCall = (Functions + Lpar - (ExpressionList) + Rpar).setName('Function call')

    Atom = (
        (FunctionCall | Field | Genotype | Number | String) | (Lpar + Expression + Rpar)
    )
    Factor = Forward()
    Factor <<= Atom + (Factor)[...]
    Term = Factor + (ArithOp + Factor)[...]
    Expression <<= Term[...]

    NullColumnCheck = ((VCFColumn | OtherColumns) + IsNull).setName('Null column check')
    ExpressionComparison = (FixUp | (Expression + Comparison + Expression)).setName('Comparison')

    # Order is important here: most to least specific
    Condition = (
        NullColumnCheck
        | ExpressionComparison
        | Field
        | Expression
    )

    return (Condition + ZeroOrMore(Logical + Condition))

def generate_missingness_parser():
    ''' Missingness parser definition '''
    Number = pyparsing_common.number
    Number.setName('Number')

    Comparison = Literal("==").setParseAction(replaceWith('=')) | oneOf("= > >= <= < !=")
    Comparison.setName('Comparison')

    VariantMiss = (
        CaselessKeyword('F_MISS') |
        CaselessKeyword('F_MISSING').setParseAction(replaceWith('F_MISS')) |
        CaselessKeyword('N_MISS') |
        CaselessKeyword('N_MISSING').setParseAction(replaceWith('N_MISS'))
    )
    VariantMiss.setParseAction(lambda s, l, t: [f"""CAST(infos['{t[0]}'] AS DOUBLE)"""])

    SampleMiss = (
        CaselessKeyword('S_MISS') |
        CaselessKeyword('S_MISSING').setParseAction(replaceWith('S_MISS')) |
        CaselessKeyword('SAMPLE_MISS').setParseAction(replaceWith('S_MISS')) |
        CaselessKeyword('SAMPLE_MISSING').setParseAction(replaceWith('S_MISS'))
    )

    Missing = (VariantMiss | SampleMiss)
    Missing.setName('F_MISS | N_MISS | S_MISS')

    return (Missing + Comparison + Number)

# Global parser instances are generated once

VDB_PARSER = generate_vdb_parser()
MISSINGNESS_PARSER = generate_missingness_parser()

def parser(query, parser_type='vdb'):
    ''' Parse a user's query string and return the equivalent SQL query suffix '''
    try:
        if parser_type == 'vdb':
            result = VDB_PARSER.parseString(query, parseAll=True)
        elif parser_type == 'missingness':
            result = MISSINGNESS_PARSER.parseString(query, parseAll=True)
        else:
            raise SystemExit(f"Invalid parser_type: {parser_type}")

        return ' '.join([str(f) for f in result])

    except ParseException as err:
        err.syntax = f'{err.line}\n{" " * (err.column - 1) + "^"}'
        raise err
