import sqlite3


class GeneDb(object):

    def __init__(self, create=False):
        self.connection = sqlite3.connect("genes.db")
        self.cursor = connection.cursor()

    def create_table(self):
        sql_command = """CREATE TABLE genes (
                id INTEGER PRIMARY KEY,
                bin INTEGER,
                name VARCHAR(20),
                name2 VARCHAR(20),
                chrom VARCHAR(5),
                strand CHAR(1),
                start_pos INTEGER,
                end_pos INTEGER);"""

        self.cursor.execute(sql_command)

    def populate_table(self, import_file):
        with open(import_file, 'r') as of:
            for line in of:
                if line[0] == '#':
                    continue
                item = line.split()
                format_str = ('INSERT INTO genes (id, bin, name, name2, chrom, strand, start_pos, end_pos) '
                              'VALUES (NULL, "{bin}", "{name}", "{name2}", "{chrom}", "{strand}", "{start}", "{end}");')
                sql_command = format_str.format(bin=item[0], name=item[1], name2=item[12], chrom=item[2],
                                                strand=item[3], start=item[4], end=item[5])
                self.cursor.execute(sql_command)
        return True

    def select_from_table(self, chromosome, position):
        # nope
        chromosome = "chr" + chromosome
        select_str = ("SELECT name2 FROM genes WHERE chrom = \'{target_chromosome}\' "
                      "AND start_pos <= {target_position} AND end_pos >= {target_position};")
        sql_command = select_str.format(target_chromosome=chromosome, target_position=position)
        cursor.execute(sql_command)
        result = cursor.fetchall()
        mystring = ""
        for r in result:
            mystring += (str(r)[3:-3] + " ")
        return mystring[:-1]


def test():
    # not how we should test
    print select_from_table("1", 75203726)
