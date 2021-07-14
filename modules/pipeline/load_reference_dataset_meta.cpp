#include "modules/pipeline/load_reference_dataset_meta.h"

#include "modules/io/log.h"
#include "modules/pipeline/dataset_meta.h" 
#include "modules/pipeline/dataset_path.h" 

void load_reference_dataset_meta()
{
	path ref_root(CONF_S(path_reference_base));
	if (!ref_root.exists()) {
		ref_root.mkdir();
	}

	dataset_path ref_root_path("/api/reference");

	auto ref_root_meta = ref_root_path.meta();
	if (!ref_root_meta.exists()) {
		ref_root_meta.mkdir();
	}

	path share_reference(CONF_S(reference_path));
	auto sources = share_reference.list();
	for (const auto& source : sources)
	{
		try
		{
			dataset_meta dm;
			auto ref = share_reference.append(source).append("dataset_meta");
			ref.json_get(dm);
			dataset_path rp(printstring("/api/reference/%s", source.c_str()));	
			rp.remove();
			rp.create(dm);
		}
		catch (const io_exception& e)
		{
			SPLOG("%s", e.message().c_str());
		}
	}
}

