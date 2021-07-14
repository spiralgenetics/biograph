#pragma once

#include "modules/bio_format/exporter.h"

class histogram_exporter : public exporter
{
public:
	histogram_exporter(writable& byte_sink);

	void write(const std::string& key, const std::string& value) override;
};
