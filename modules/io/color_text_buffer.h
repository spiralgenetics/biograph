
#include "modules/io/io.h"
#include <map>

namespace color {
	const uint32_t red = 0xff0000;
	const uint32_t green = 0x00ff00; 
	const uint32_t blue = 0x0000ff;
	const uint32_t white = 0xffffff; 
	const uint32_t black = 0x000000;
	const uint32_t grey = 0x808080; 
}

class color_text_buffer {
public:
	// Construct an empty text buffer
	color_text_buffer();
	// Set's position (negative values are allowed)
	void set_position(int x, int y);
	// Set's color
	void set_color(uint32_t rgb);
	// Adds text
	void print(const char* fmt, ...) __attribute__((format(printf,2,3)));
	// Write it as html
	void render_as_html(writable& out);
private:
	struct char_t { uint32_t color; char c; };
	typedef std::map<int, char_t> row_t;
	typedef std::map<int, row_t> screen_t;
	uint32_t m_color;
	int m_xpos;
	int m_ypos;
	screen_t m_screen;
};

