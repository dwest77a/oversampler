import json
INPUTS = {
	'iasi':[
		'/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/iasi_l2/valid_l2a/v00.00/{}/{}/{}',
		'/gws/pw/j05/rsg_share/public/projects/ims/data/lv2/output_ghg_cv9_v123_lam_nat_nfo33_fgsnwp_nobc_newbc_rbc_ram4_rem20-ncam_rcl5_raer6st_rmg_oap_cnv2_ixam7-8_o3f110_swir2/metopa/{}/{}/{}'
	],
	'cris':[
		'/gws/pw/j05/rsg_share/public/transfer/barry/6ad76b04-b3df-11eb-a4d3-024ad8e814ad/cris_l2/valid_l2a/v00.00c/{}/{}/{}',
		''
	],
    'output':[
        '',
        ''
    ]
}

filename = '../base_validation/file_config.json'

f = open(filename,"w")
f.write(json.dumps(INPUTS))
f.close()