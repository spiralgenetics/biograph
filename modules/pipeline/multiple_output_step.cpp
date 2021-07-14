#include "modules/pipeline/multiple_output_step.h"
#include "modules/io/log.h"

void multiple_output_step::finalize(
	dataset_meta& out,
	const std::string& output,
	const std::vector<dataset_meta>& inputs,
	const std::string& options) const
{
	SPLOG("multiple_output_step::finalize> #%zu", out.get_output_index());
	std::vector<manifest> output_manifests;
	json_deserialize(output_manifests, output);
	out.the_manifest = output_manifests.at(out.get_output_index());
	update_metadata(out, inputs, options);
}
