#pragma once

#include <string>
#include "modules/mapred/manifest.h"

class scoped_temp_file
{
public:
	scoped_temp_file();
	explicit scoped_temp_file(const std::string& tmpl);
	~scoped_temp_file();

	std::string path() const
	{
		return m_path;
	}

private:
	std::string m_path;
        int m_fd;
};

// for pipe_mapper
struct temp_file_spec
{
        TRANSFER_OBJECT {
                VERSION(0);

                FIELD(arg_index, TF_STRICT);
                FIELD(data, TF_STRICT);
                FIELD(exporter_type, TF_STRICT);
                FIELD(ex_im_porter_data);
        }

        temp_file_spec() : arg_index(std::numeric_limits<unsigned short>::max()) {}
        temp_file_spec(unsigned short an_arg_index, const manifest& the_data,
                                        const std::string& the_exporter_type, const std::string& the_ex_im_porter_data)
                : arg_index(an_arg_index), data(the_data), exporter_type(the_exporter_type),
                ex_im_porter_data(the_ex_im_porter_data) {}

        unsigned short arg_index;
        manifest data;
        std::string exporter_type;
        std::string ex_im_porter_data;
};
