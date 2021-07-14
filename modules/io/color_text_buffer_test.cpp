#include "modules/test/test_utils.h"
#include "modules/io/color_text_buffer.h"
#include "modules/io/file_io.h"
#include <gtest/gtest.h>

TEST(color_text_buffer, one) 
{
	color_text_buffer ctb;
	ctb.set_position(0, 3);
	ctb.print("Row 3");
	ctb.set_position(-1, 2);
	ctb.print("Row 2, minus 1");
	ctb.set_position(5, 5);
	ctb.set_color(color::red);
	ctb.print("Red");
	ctb.set_position(6, 6);
	ctb.set_color(color::green);
	ctb.print("Green");
	ctb.set_position(7, 7);
	ctb.set_color(color::blue);
	ctb.print("Blue");
	ctb.set_position(0, 0);
	ctb.print("This will be overwritten");
	ctb.set_position(13, 0);
	ctb.set_color(color::red);
	ctb.print("Awesome!   ");
	file_writer out(make_path("color_text_buffer.html"));
	out.print("<html><body>\n");
	ctb.render_as_html(out);
	out.print("</body></html>\n");
}

