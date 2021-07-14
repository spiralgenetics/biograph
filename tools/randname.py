#!/usr/bin/env python3
#
# Random improbable name generator
#
# Syllables borrowed from https://github.com/mouse-reeve/fruit
#

import random
import string
import sys

left = (
    'Ace',
    'Alm',
    'Ap',
    'Ate',
    'Avo',
    'Ba',
    'Bil',
    'Bir',
    'Cai',
    'Can',
    'Car',
    'Cha',
    'Cher',
    'Cit',
    'Cur',
    'Fi',
    'Go',
    'Gra',
    'Gua',
    'Ja',
    'Jac',
    'Jam',
    'Ju',
    'Ka',
    'Ki',
    'Le',
    'Li',
    'Lin',
    'Lo',
    'Lon',
    'Ly',
    'Ma',
    'Ma',
    'Mar',
    'Med',
    'Mel',
    'Mor',
    'Pa',
    'Paw',
    'Pe',
    'Pep',
    'Per',
    'Plu',
    'Po',
    'Pom',
    'Pru',
    'Qui',
    'Sap',
    'Tam',
    'Tan',
    'Tang',
    'To',
)

middle = (
    'amb',
    'ang',
    'ar',
    'ba',
    'cad',
    'co',
    'ero',
    'er',
    'ger',
    'go',
    'gran',
    'ist',
    'ju',
    'ma',
    'mi',
    'mo',
    'moy',
    'na',
    'pa',
    'po',
    'ri',
    'rim',
    'ru',
    'sim',
    'tel',
    'yo',
)

right = (
    'ach',
    'ang',
    'ar',
    'at',
    'ate',
    'be',
    'bul',
    'ce',
    'chee',
    'cot',
    'go',
    'ngo',
    'ig',
    'inda',
    'ine',
    'ito',
    'kin',
    'la',
    'lon',
    'mon',
    'ne',
    'nge',
    'ola',
    'on',
    'ot',
    'ond',
    'osu',
    'oupe',
    'out',
    'paw',
    'pe',
    'per',
    'ple',
    'pok',
    'quat',
    'rant',
    'rind',
    'rine',
    'ron',
    'ry',
    'te',
    'to',
    'tron',
    'ula',
    'wi',
    'ya',
)

def get_word():
    return random.choice(left) + ''.join(random.sample(middle, random.randint(0,1))) + random.choice(right)

def get_random_name(seed=None):
    if seed:
        random.seed(seed)

    return f"{get_word()} {random.choice(string.ascii_uppercase)}. {get_word()}"

if __name__ == '__main__':
    try:
        print(get_random_name(sys.argv[1]))
    except IndexError:
        print(get_random_name())
