#pragma once

#include "modules/pipeline/steptype.h"

// Intermediate base class for steps with more than one output.  Inherit from this class
// to get the finalize that deals with the output being a vector of manifests instead of
// a single manifest.
class multiple_output_step : public steptype
{
public:
	void finalize(
		dataset_meta& out, 
		const std::string& output, 
		const std::vector<dataset_meta>& inputs, 
		const std::string& options) const override;
};
