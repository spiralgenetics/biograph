#include "modules/pipeline/steptype.h"
#include "modules/pipeline/dataset_path.h" 

void validate_one_input(const std::string& input, const steptype::param& input_type, const std::string& stepname)
{
	dataset_path ds(input);
	direntry de = ds.stat();
	if (de.directory) {
		throw io_exception(printstring("Cannot run command on a directory: %s", ds.friendly().c_str()));
	}

	if (de.type != input_type.type) {
		throw io_exception(printstring("Wrong argument type for step of type '%s', type of file %s was %s, should be %s",
			stepname.c_str(),
			ds.friendly().c_str(),
			de.type->id.c_str(),
			input_type.type->id.c_str()
		));
	}
	
	dataset_meta input_dataset_meta;
	ds.load(input_dataset_meta);
	if (input_dataset_meta.in_progress) {
		throw io_exception(
			boost::format("Cannot use input %1% because it is still being generated.") % ds.friendly()
		);
	}
}

void strict_validation(const steptype& steptype, const std::vector<std::string>& inputs)
{
	if (inputs.size() != steptype.inputs.size()) {
		throw io_exception(printstring("Incorrect number of inputs for step of type '%s', should be %ld, was %ld", 
			steptype.name.c_str(), 
			steptype.inputs.size(),
			inputs.size()
		));
	}

	for (size_t i = 0; i < inputs.size(); i++) {
		try {
			validate_one_input(inputs[i], steptype.inputs[i], steptype.name);
		}
		catch (const io_exception& io) {
			throw io_exception(printstring("Input %ld: %s", i + 1, io.message().c_str()));
		}
	}
}

void many_same_type_validation(const steptype& steptype, const std::vector<std::string>& inputs)
{
	if (inputs.empty()) {
		throw io_exception(printstring("Missing inputs for step of type '%s'", 
			steptype.name.c_str()
		));
	}

	if (steptype.inputs.size() != 1) {
		throw io_exception(printstring("Invalid Steptype Definition: Incorrect number of inputs for type '%s'", 
			steptype.name.c_str()
		));
	}

	steptype::param input_type = steptype.inputs[0];
	for (size_t i = 0; i < inputs.size(); i++) {
		try {
			validate_one_input(inputs[i], input_type, steptype.name);
		}
		catch (const io_exception& io) {
			throw io_exception(printstring("Input %lu: %s: Error: %s", i + 1, 
				inputs[i].c_str(), 
				io.message().c_str()
			));
		}
	}
}

void no_check_validation(const steptype& steptype, const std::vector<std::string>& inputs) {}

std::map<int, std::function<void (const steptype&, const std::vector<std::string>& )> > g_validation_functions
{
	{ steptype::STRICT, strict_validation },
	{ steptype::MANY_SAME_TYPE, many_same_type_validation },
	{ steptype::NO_CHECK, no_check_validation }
};

void steptype::validate_input(const std::vector<std::string>& inputs) const
{
	g_validation_functions[input_validation_policy](*this, inputs);
}

void steptype::add_input(
	const char* type, 
	const char* name, 
	const std::vector<std::string>& sort_keys,
	bool is_optional)
{
	param p;
	p.type = datatype_registry::find(type);
	p.name = name;
	p.sort_keys = sort_keys;
	p.is_optional = is_optional;
	inputs.push_back(p);
}

void steptype::add_output(
	const char* type, 
	const char* name, 
	const std::vector<std::string>& sort_keys,
	bool is_optional)
{
	param p;
	p.type = datatype_registry::find(type);
	p.name = name;
	p.sort_keys = sort_keys;
	p.is_optional = is_optional;
	outputs.push_back(p);
}

void steptype::update_metadata(
	dataset_meta& out, 
	const std::vector<dataset_meta>& inputs, 
	const std::string& options) const
{
	for (const auto& input : inputs) {
		out.the_manifest.merge_tags(input.the_manifest);
	}
	update_step_metadata(out, options);
}

void steptype::finalize(
	dataset_meta& out, 
	const std::string& output, 
	const std::vector<dataset_meta>& inputs, 
	const std::string& options) const
{
	SPLOG("steptype::finalize> %s", output.c_str());
	json_deserialize(out.the_manifest, output);
	update_metadata(out, inputs, options);
}

void steptype::set_runtime_in_metadata(dataset_meta& out, time_t start_time)
{
	out.the_manifest.metadata().set_runtime(start_time);
}

void steptype::update_step_metadata(dataset_meta& out, const std::string& options) const
{
	out.the_manifest.metadata().set_creation_time_now();
	out.the_manifest.metadata().set_options(id, options);
}
