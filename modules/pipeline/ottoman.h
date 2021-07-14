#pragma once

#include <memory>

std::string ottoman_url();

class ottoman_impl;

class ottoman_server
{
public:
	ottoman_server();
	~ottoman_server();

private:
	std::unique_ptr<ottoman_impl> m_impl;
};
