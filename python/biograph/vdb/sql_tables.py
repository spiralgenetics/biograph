"""

Athena SQL table definitions.

NOTE: make sure that any changes here are reflected in the parquet schemas defined
in parquet_variants.py and parquet_anno.py.

"""

def create_table_sql(table_type, table_name=None):
    """
    Generates SQL DDL to create a variants table with the given name.
    Returns a query which has a single named parameter, 'location', which should
    refer to the s3 location of the table data.
    """

    # The table name is usually the same as its type, but this lets us
    # easily create new tables of a given type for development
    if not table_name:
        table_name = table_type

    if "`" in table_name:
        raise RuntimeError(f"Table name {table_name} may not contain backquotes")
    try:
        return {
            "variants":
                f"""
                    CREATE EXTERNAL TABLE IF NOT EXISTS `{table_name}` (
                        `chrom` STRING,
                        `spans` ARRAY<BIGINT>,
                        `reflen` BIGINT,
                        `pos` BIGINT,
                        `varend` BIGINT,
                        `varid` STRING,
                        `ref` STRING,
                        `alt` STRING,
                        `qual` FLOAT,
                        `filt` STRING,
                        `info` MAP<STRING, STRING>,
                        `sample` MAP<STRING, STRING>
                    )
                    PARTITIONED BY (
                        `sample_name` STRING,
                        `build` STRING,
                        `aid` STRING
                    )
                    STORED AS PARQUET
                    LOCATION %(location)s
                    ;
                """,

            "headers":
                f"""
                    CREATE EXTERNAL TABLE IF NOT EXISTS `{table_name}` (
                        `refname` STRING,
                        `refhash` STRING,
                        `header` STRING,
                        `description` STRING,
                        `imported_on` TIMESTAMP,
                        `variant_count` BIGINT,
                        `filename` STRING
                    )
                    PARTITIONED BY (
                        `sample_name` STRING,
                        `build` STRING,
                        `aid` STRING
                    )
                    ROW FORMAT SERDE 'org.openx.data.jsonserde.JsonSerDe'
                    LOCATION %(location)s
                    ;
                """,

            "study_variants":
                f"""
                    CREATE EXTERNAL TABLE IF NOT EXISTS `{table_name}` (
                        `chrom` STRING,
                        `spans` ARRAY<BIGINT>,
                        `reflen` BIGINT,
                        `pos` BIGINT,
                        `varend` BIGINT,
                        `varid` STRING,
                        `ref` STRING,
                        `alt` STRING,
                        `qual` FLOAT,
                        `filt` STRING,
                        `info` MAP<STRING, STRING>,
                        `sample` MAP<STRING, STRING>
                    )
                    PARTITIONED BY (
                        `study_name` STRING,
                        `checkpoint` BIGINT,
                        `sample_name` STRING,
                        `aid` STRING
                    )
                    STORED AS PARQUET
                    LOCATION %(location)s
                    ;
                """,

            "study_merged":
                f"""
                    CREATE EXTERNAL TABLE IF NOT EXISTS `{table_name}` (
                        `chrom` STRING,
                        `spans` ARRAY<BIGINT>,
                        `reflen` BIGINT,
                        `pos` BIGINT,
                        `varend` BIGINT,
                        `varid` STRING,
                        `ref` STRING,
                        `alt` STRING,
                        `qual` FLOAT,
                        `filters` ARRAY<STRING>,
                        `infos` MAP<STRING, STRING>,
                        `samples` MAP<STRING, MAP<STRING, STRING>>
                    )
                    PARTITIONED BY (
                        `study_name` STRING,
                        `checkpoint` BIGINT
                    )
                    STORED AS PARQUET
                    LOCATION %(location)s
                    ;
                """,

            "study_meta":
                f"""
                    CREATE EXTERNAL TABLE IF NOT EXISTS `{table_name}` (
                        `key` STRING,
                        `value` STRING,
                        `dtype` STRING
                    )
                    PARTITIONED BY (
                        `study_name` STRING
                    )
                    ROW FORMAT DELIMITED
                    FIELDS TERMINATED BY '\\t'
                    ESCAPED BY '\\\\'
                    LINES TERMINATED BY '\\n'
                    LOCATION %(location)s
                    ;
                """,

            "anno":
                f"""
                    CREATE EXTERNAL TABLE IF NOT EXISTS `{table_name}` (
                        `chrom` STRING,
                        `spans` ARRAY<BIGINT>,
                        `reflen` BIGINT,
                        `pos` BIGINT,
                        `varend` BIGINT,
                        `varid` STRING,
                        `ref` STRING,
                        `alt` STRING,
                        `qual` FLOAT,
                        `filt` STRING,
                        `info` MAP<STRING, STRING>,
                        `source` STRING,
                        `feature` STRING,
                        `score` FLOAT,
                        `frame` STRING,
                        `strand` STRING,
                        `attributes` MAP<STRING, STRING>
                    )
                    PARTITIONED BY (
                        `build` STRING,
                        `anno_name` STRING,
                        `version` STRING,
                        `aid` STRING
                    )
                    STORED AS PARQUET
                    LOCATION %(location)s
                    ;
                """,

            "anno_meta":
                f"""
                    CREATE EXTERNAL TABLE IF NOT EXISTS `{table_name}` (
                        `refname` STRING,
                        `refhash` STRING,
                        `header` STRING,
                        `description` STRING,
                        `imported_on` TIMESTAMP,
                        `variant_count` BIGINT,
                        `filename` STRING
                    )
                    PARTITIONED BY (
                        `build` STRING,
                        `anno_name` STRING,
                        `version` STRING,
                        `aid` STRING
                    )
                    ROW FORMAT SERDE 'org.openx.data.jsonserde.JsonSerDe'
                    LOCATION %(location)s
                    ;
                """,

        }[table_type]

    except KeyError:
        raise SystemExit(f"Cannot create table of unknown type: {table_type}")
