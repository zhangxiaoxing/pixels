from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import pandas as pd
from allensdk.api.queries.ontologies_api import OntologiesApi

# The manifest file is a simple JSON file that keeps track of all of
# the data that has already been downloaded onto the hard drives.
# If you supply a relative path, it is assumed to be relative to your
# current working directory.

def get_conn_target():
    mcc = MouseConnectivityCache(manifest_file='connectivity/mouse_connectivity_manifest.json')
    structure_tree = mcc.get_structure_tree()
    conn_target = structure_tree.get_structures_by_set_id([184527634])
    return [x['id'] for x in conn_target]

def get_conn_inj():
    mcc = MouseConnectivityCache(manifest_file='connectivity/mouse_connectivity_manifest.json')
    structure_tree = mcc.get_structure_tree()
    conn_target = structure_tree.get_structures_by_set_id([114512891])
    return [x['id'] for x in conn_target]


if __name__=='__main__':
    mcc = MouseConnectivityCache(manifest_file='connectivity/mouse_connectivity_manifest.json')
    structure_tree = mcc.get_structure_tree()
    oapi = OntologiesApi()

    # get the ids of all the structure sets in the tree
    structure_set_ids = structure_tree.get_structure_sets()

    # query the API for information on those structure sets
    t=pd.DataFrame(oapi.get_structure_sets(structure_set_ids))

    # On the connectivity atlas web site, you'll see that we show most of our data at a fairly coarse structure level. We did this by creating a structure set of ~300 structures, which we call the "summary structures". We can use the structure tree to get all of the structures in this set:

    # From the above table, "Mouse Connectivity - Summary" has id 167587189
    summary_structures = structure_tree.get_structures_by_set_id([3])
    t=pd.DataFrame(summary_structures)