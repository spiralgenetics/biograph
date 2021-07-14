#include "modules/io/pass_thru.h"
#include <stdarg.h>

void pass_thru_writable::print(const char* format, ...)
{
	va_list args;
	va_start(args, format);
	m_in.print(format, args);
	va_end(args);
}
