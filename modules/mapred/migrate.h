#pragma once

#include "modules/mapred/manifest.h" 
#include <vector> 

typedef std::function<void(manifest& dataset)> migration_f;

void migrate_000_001(manifest& dataset);
void migrate_001_002(manifest& dataset);
void migrate_002_003(manifest& dataset);
void migrate_003_004(manifest& dataset);

const std::vector<migration_f> migrations {
	migrate_000_001, 
	migrate_001_002,
	migrate_002_003,
	migrate_003_004
};
