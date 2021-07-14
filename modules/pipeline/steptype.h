#pragma once

#include "modules/pipeline/restful_registry.h"
#include "modules/pipeline/datatype.h"
#include "modules/pipeline/dataset_meta.h"
#include "modules/mapred/task.h"
#include "modules/io/transfer_object.h"

class steptype
{
public:
	struct param
	{
		TRANSFER_OBJECT
		{
			VERSION(0);
			FIELD(type, TF_STRICT); 
			FIELD(name, TF_STRICT); 
			FIELD(sort_keys, TF_STRICT);
			FIELD(is_optional, false);
		}

		datatype_ref type;
		std::string name;
		std::vector<std::string> sort_keys;
		bool is_optional;
	};
	
	TRANSFER_OBJECT
	{
		VERSION(0);
		FIELD(id, TF_STRICT); 
		FIELD(url, TF_STRICT); 
		FIELD(name, TF_STRICT); 
		FIELD(description, TF_STRICT); 
		FIELD(input_validation_policy, TF_STRICT); 
		FIELD(output_validation_policy, TF_STRICT); 
		FIELD(inputs, TF_STRICT); 
		FIELD(outputs, TF_STRICT); 
	}

	std::string id;
	std::string url;
	std::string name;
	std::string description;
	std::vector<param> inputs;
	std::vector<param> outputs;
	int input_validation_policy = ILLEGAL_INPUT_VALIDATION;
	int output_validation_policy = ILLEGAL_OUTPUT_VALIDATION;

	enum input_validation_policy_e 
	{
		ILLEGAL_INPUT_VALIDATION = 0, // A default policy that will throw if not overridden.
		STRICT, // task inputs must match in number and types with the ones in this steptype
		MANY_SAME_TYPE, // task inputs can vary in number (>0) but must all be of same type (defined in inputs[0] )
		NO_CHECK // does not perform any checking
	};

	enum output_validation_policy_e 
	{
		ILLEGAL_OUTPUT_VALIDATION = 0, // A default policy that will throw if not overridden.
		PIPELINE_STEP, // must have only one output. Must be a file that does not exist yet
		NO_OUTPUT // step produces no output
	};
	
	void add_input(
		const char* type, 
		const char* name, 
		const std::vector<std::string>& sort_keys = std::vector<std::string>(),
		bool is_optional = false
	);

	void add_output(
		const char* type, 
		const char* name, 
		const std::vector<std::string>& sort_keys = std::vector<std::string>(),
		bool is_optional = false
	);
	
	unsigned get_optional_output_count() const
	{
		return std::count_if(outputs.cbegin(), outputs.cend(),
			[](const param& output_param) { return output_param.is_optional; }
		);
	}

	void validate_input(const std::vector<std::string>& inputs) const;

	// Verify both options and inputs (throws if failure), and compute goal machines
	virtual void validate(
		const std::vector<std::string>& inputs, 
		const std::string& options
	) const = 0;

	virtual std::unique_ptr<task> create_task(
		const std::vector<dataset_meta>& inputs, 
		const std::string& options
	) const = 0;

	// The finalization is separated into two parts: update_metadata transfers metadata from step input to outputs and is called from
	// the pre-run validation (to check inputs and outputs) and from the finalize method to update the actual output metadata.  The
	// finalize method calls update_metadata and then writes the actual output manifests.
	//
	// Note that the base class steptype::update metadata updates several metadata items including the
	// time of the run and the software version that ran.  If you override update_metadata, be sure to
	// either call the base class or implement any of the base class functionality that you desire (if
	// any) in the derived class.
	virtual void update_metadata(
		dataset_meta& out, 
		const std::vector<dataset_meta>& inputs, 
		const std::string& options
	) const;
	
	virtual void finalize(
		dataset_meta& out, 
		const std::string& output, 
		const std::vector<dataset_meta>& inputs, 
		const std::string& options
	) const;
	
	void set_runtime_in_metadata(dataset_meta& out, time_t start_time);
	void update_step_metadata(dataset_meta& outconst, const std::string& options) const;
};

typedef restful_registry<steptype> steptype_registry;
typedef steptype_registry::ref_type steptype_ref;
