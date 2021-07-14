
#include "modules/io/color_text_buffer.h"
#include "modules/io/log.h"
#include <stdarg.h>
#include <stdio.h>
#include <limits>

color_text_buffer::color_text_buffer() 
	: m_color(color::white) 
	, m_xpos(0)
	, m_ypos(0)
{}

void color_text_buffer::set_position(int x, int y) 
{
	//SPLOG("Position: %d, %d", x, y);
	m_xpos = x;
	m_ypos = y;
}

void color_text_buffer::set_color(uint32_t rgb) 
{
	//SPLOG("Color: %06x", rgb);
	m_color = rgb;
}

void color_text_buffer::print(const char* fmt, ...)
{
        va_list vl;
        va_start(vl, fmt);

        char* buf;
        int len = vasprintf(&buf, fmt, vl);
	//SPLOG("%s", buf);

	char_t c;
	c.color = m_color;
	for(size_t i = 0; i < (size_t) len; i++)
	{
		c.c = buf[i];
		m_screen[m_ypos][m_xpos++] = c;
	}

        free(buf);
        va_end(vl);
}

void color_text_buffer::render_as_html(writable& out) 
{
	int min_left = std::numeric_limits<int>::max();
	for(const auto& kvp : m_screen) {
		min_left = std::min(min_left, kvp.second.begin()->first);
	}
	uint32_t cur_color = color::white;
	out.print("<pre><span style=\"background-color: #%6x\">   ", cur_color);
	int cur_y = m_screen.begin()->first;
	for(const auto& row : m_screen) {
		while (cur_y < row.first) { out.print("   \n   "); cur_y++; }
		int cur_x = min_left;
		for(const auto& col : row.second) {
			while(cur_x < col.first) { out.print(" "); cur_x++; }
			if (cur_color != col.second.color) {
				cur_color = col.second.color;
				out.print("</span><span style=\"background-color: #%06x\">", cur_color);
			}
			out.print("%c", col.second.c);
			cur_x++;
		}
	}
	out.print("   \n</span></pre>");
}



