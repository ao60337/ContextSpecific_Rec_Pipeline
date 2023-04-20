
import pandas as pd
from Tissue_specific_Reconstruction_Pipeline import DATA_PATH

media = pd.read_excel(rf"{DATA_PATH}/media.xlsx", index_col=0, sheet_name=None)['media_with_starch'].to_dict(orient='index')
MEDIUM_CONDITIONS = {'f2 medium': {key: (value['LB'], value['UB']) for key, value in media.items()}}

UPTAKE_DRAINS = {'f2 medium': {"EX_C00014__dra", "EX_C00059__dra", "EX_C01330__dra", "EX_C00011__dra", "EX_C14818__dra", "EX_C00080__dra", "EX_C00001__dra", "EX_C00305__dra", "EX_C01327__dra", "EX_C00244__dra", "EX_C00009__dra",
                                                                       "EX_C00007__dra", "EX_C00205__dra", "EX_C00378__dra", "EX_C00120__dra", "EX_C02823__dra", 'DM_C00369__chlo'
                                                                       }}

MEDIUM_METABOLITES = {'f2 medium': [e.split("__")[0].replace("EX_", "") for e in MEDIUM_CONDITIONS['f2 medium']]}
