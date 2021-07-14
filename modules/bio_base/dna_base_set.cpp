
#include "modules/bio_base/dna_base_set.h"
#include "modules/io/utils.h"

#include <map>
#include <algorithm>
#include <boost/bind.hpp>

// Private loookup table class, uses static construction to fill table
// Implements IUPAC Nucleic Acid Code
// As describe here: http://www.dna.affrc.go.jp/misc/MPsrch/InfoIUPAC.html
// Note: we ignore the 'U' case, since this is for DNA

class code_table
{
public:
	code_table()
	{
		dna_base_set a = dna_base('A');
		dna_base_set c = dna_base('C');
		dna_base_set g = dna_base('G');
		dna_base_set t = dna_base('T');
		add('-', dna_base_set());
		add('A', a);
		add('C', c);
		add('G', g);
		add('T', t);
		add('M', a | c);
		add('R', a | g);
		add('W', a | t);
		add('S', c | g);
		add('Y', c | t);
		add('K', g | t);
		add('V', a | c | g);
		add('H', a | c | t);
		add('D', a | g | t);
		add('B', c | g | t);
		add('N', a | c | g | t);
	}
	char lookup(const dna_base_set& b) const { return m_reverse.find(b)->second; }
	const dna_base_set& lookup(char c) const { return m_forward.find(c)->second; }
private:
	void add(char c, dna_base_set b) 
	{ 
		m_forward.insert(std::make_pair(c,b));
		m_reverse.insert(std::make_pair(b,c));
	}
	std::map<char, dna_base_set> m_forward;
	std::map<dna_base_set, char> m_reverse;
};

static code_table the_table;  // Global for static construction

dna_base_set::dna_base_set(char c)
	: m_set(the_table.lookup(c).m_set)
{}

char dna_base_set::as_code()
{
	return the_table.lookup(*this);
}

void print_a_base(const dna_base& base, std::string& result, const std::string& separator)
{
	result = printstring("%s%c%s", result.c_str(), (char)base, separator.c_str());
}

std::string dna_base_set::as_list(char sep)
{
	std::string str_sep;
	if (sep) str_sep.push_back(sep);
	std::string result;
	dna_base_set::reduce(boost::bind(print_a_base, _1, _2, str_sep), *this, result);
	if (sep && result.size()) return result.substr(0, result.size() - 1);
	return result;
}

//                                             ABCD..GH..K.MN...RST.VW.Y.
static constexpr char k_iupac_complements[] = "TVGH..CD..M.KN...YSA.BW.R.";
static constexpr int k_num_iupac_complements =
    sizeof(k_iupac_complements) - 1 /* don't count the \0 at end */;

void reverse_complement_iupac_string(std::string& iupac_string) {
  std::reverse(iupac_string.begin(), iupac_string.end());
  std::transform(iupac_string.begin(), iupac_string.end(), iupac_string.begin(),
                 [](const char& c) -> char {
                   int idx = int(c) - 'A';
                   if (idx < 0 || idx >= k_num_iupac_complements) {
                     return c;
                   } else {
                     return k_iupac_complements[idx];
                   }
                 });
}
