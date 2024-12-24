#include "search.h"
#include "settings.h"
#include "tt.h"
#include "types.h"

struct settings settings, delayedSettings;

SMALL void process_delayed_settings(void) {
    if (delayedSettings.ttSize != settings.ttSize)
    {
        tt_free();
        settings.ttSize = delayedSettings.ttSize;
        tt_allocate(settings.ttSize);
    }

    if (delayedSettings.clear)
    {
        delayedSettings.clear = false;
        search_clear();
    }
}
