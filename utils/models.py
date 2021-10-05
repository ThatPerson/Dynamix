mds = {
	'smf': {
		'n': 2,
		'p': ['tau', 'S2']
	}, 
	'smft': {
		'n': 3,
		'p': ['tau', 'S2', 'Ea']
	}, 
	'emf': {
		'n': 3,
		'p': ['taus', 'S2s', 'tauf']
	}, 
	'emft': {
		'n': 5,
		'p': ['taus', 'S2s', 'tauf', 'Eas', 'Eaf']
	}, 
	'demf': {
		'n': 4,
		'p': ['taus', 'S2s', 'tauf', 'S2f']
	}, 
	'demft': {
		'n': 6,
		'p': ['taus', 'S2s', 'tauf', 'S2f', 'Eas', 'Eaf']
	}, 
	'rdemft': {
		'n': 8,
		'p': ['taus', 'S2s', 'tauf', 'S2f', 'Eas', 'Eaf', 'papbS', 'kex']
	}, 
	'sdemft': {
		'n': 8,
		'p': ['taus', 'S2s', 'tauf', 'S2f', 'Eas', 'Eaf', 'dS2s', 'dS2f']
	}, 
	'sdemf': {
		'n': 6,
		'p': ['taus', 'S2s', 'tauf', 'S2f', 'dS2s', 'dS2f']
	}, 
	'gdemf': {
		'n': 6,
		'p': ['taus', 'S2s', 'tauf', 'S2f', 'GC0', 'Gtr']
	},
	'gsdemf': {
		'n': 6,
		'p': ['taus', 'S2s', 'tauf', 'S2f', 'GC0', 'Gtr']
	},
	'gsmf': {
		'n': 4,
		'p': ['tau', 'S2', 'GC0', 'Gtr']
	},
	'demft_no600': {
		'n': 6,
		'p': ['taus', 'S2s', 'tauf', 'S2f', 'Eas', 'Eaf']
	}, 
	'demft600': {
		'n': 7,
		'p': ['taus', 'S2s', 'tauf', 'S2f', 'Eas', 'Eaf', '600del']
	}, 
	'egaf': {
		'n': 6,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sGs', 'S2f']
	}, 
	'egaft': {
		'n': 8,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sGs', 'S2f', 'Eas', 'Eaf']
	}, 
	'gaf': {
		'n': 8,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sGs', 'sAf', 'sBf', 'sGf']
	}, 
	'gaft': {
		'n': 10,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sGs', 'sAf', 'sBf', 'sGf', 'Eas', 'Eaf']
	}, 
	'bgf': {
		'n': 8,
		'p': ['taus', 'tauf', 'ss', 'As', 'Bs', 'sf', 'Af', 'Bf']
	},
	'bgft': {
		'n': 10,
		'p': ['taus', 'tauf', 'ss', 'As', 'Bs', 'sf', 'Af', 'Bf', 'Eas', 'Eaf']
	},
	'aimf': {
		'n': 8,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sCs', 'sAf', 'sBf', 'sCf']
	}, 
	'aimft': {
		'n': 10,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sCs', 'sAf', 'sBf', 'sCf', 'Eas', 'Eaf']
	}, 
	'vaimf': {
		'n': 8,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sCs', 'sAf', 'sBf', 'sCf', 'alph', 'beta', 'gamm']
	}, 
	'vaimft': {
		'n': 10,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sCs', 'sAf', 'sBf', 'sCf', 'Eas', 'Eaf', 'alph', 'beta', 'gamm']
	}, 
	'vegaf': {
		'n': 9,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sGs', 'S2f', 'alph', 'beta', 'gamm']
	}, 
	'vegaft': {
		'n': 11,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sGs', 'S2f', 'Eas', 'Eaf', 'alph', 'beta', 'gamm']
	}, 
	'vgaf': {
		'n': 11,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sGs', 'sAf', 'sBf', 'sGf', 'alph', 'beta', 'gamm']
	}, 
	'vgaft': {
		'n': 13,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sGs', 'sAf', 'sBf', 'sGf', 'Eas', 'Eaf', 'alph', 'beta', 'gamm']
	}, 	
	'usmf': {
		'n': 3,
		'p': ['tau', 'S2', 'S2uf']
	}, 
	'usmft': {
		'n': 4,
		'p': ['tau', 'S2', 'Ea', 'S2uf']
	}, 
	'uemf': {
		'n': 4,
		'p': ['taus', 'S2s', 'tauf', 'S2uf']
	}, 
	'uemft': {
		'n': 6,
		'p': ['taus', 'S2s', 'tauf', 'Eas', 'Eaf', 'S2uf']
	}, 
	'udemf': {
		'n': 5,
		'p': ['taus', 'S2s', 'tauf', 'S2f', 'S2uf']
	}, 
	'udemft': {
		'n': 7,
		'p': ['taus', 'S2s', 'tauf', 'S2f', 'Eas', 'Eaf', 'S2uf']
	}, 
	'uegaf': {
		'n': 7,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sGs', 'S2f', 'S2uf']
	}, 
	'uegaft': {
		'n': 9,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sGs', 'S2f', 'Eas', 'Eaf', 'S2uf']
	}, 
	'ugaf': {
		'n': 9,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sGs', 'sAf', 'sBf', 'sGf', 'S2uf']
	}, 
	'ugaft': {
		'n': 11,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sGs', 'sAf', 'sBf', 'sGf', 'Eas', 'Eaf', 'S2uf']
	}, 
	'ubgf': {
		'n': 9,
		'p': ['taus', 'tauf', 'ss', 'As', 'Bs', 'sf', 'Af', 'Bf', 'S2uf']
	},
	'ubgft': {
		'n': 11,
		'p': ['taus', 'tauf', 'ss', 'As', 'Bs', 'sf', 'Af', 'Bf', 'Eas', 'Eaf', 'S2uf']
	},
	'uaimf': {
		'n': 9,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sCs', 'sAf', 'sBf', 'sCf', 'S2uf']
	}, 
	'uaimft': {
		'n': 11,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sCs', 'sAf', 'sBf', 'sCf', 'Eas', 'Eaf', 'S2uf']
	}, 
	'uvaimf': {
		'n': 9,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sCs', 'sAf', 'sBf', 'sCf', 'S2uf', 'alph', 'beta', 'gamm']
	}, 
	'uvaimft': {
		'n': 11,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sCs', 'sAf', 'sBf', 'sCf', 'Eas', 'Eaf', 'S2uf', 'alph', 'beta', 'gamm']
	}, 
	
	
	
	'ueaimf': {
		'n': 7,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sCs', 'S2f', 'S2uf']
	}, 
	'ueaimft': {
		'n': 9,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sCs', 'S2f', 'Eas', 'Eaf', 'S2uf']
	}, 
	'uveaimf': {
		'n': 10,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sCs', 'S2f', 'S2uf', 'alph', 'beta', 'gamm']
	}, 
	'uveaimft': {
		'n': 12,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sCs', 'S2f', 'Eas', 'Eaf', 'S2uf', 'alph', 'beta', 'gamm']
	}, 
	'eaimf': {
		'n': 6,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sCs', 'S2f']
	}, 
	'eaimft': {
		'n': 8,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sCs', 'S2f', 'Eas', 'Eaf']
	}, 
	'veaimf': {
		'n': 9,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sCs', 'S2f', 'alph', 'beta', 'gamm']
	}, 
	'veaimft': {
		'n': 11,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sCs', 'S2f', 'Eas', 'Eaf', 'alph', 'beta', 'gamm']
	}, 
	
	'baimf': {
		'n': 4,
		'p': ['tau', 'sAs', 'sBs', 'sCs']
	}, 
	'baimft': {
		'n': 6,
		'p': ['tau', 'sAs', 'sBs', 'sCs', 'Ea']
	}, 
	'vbaimf': {
		'n': 7,
		'p': ['tau', 'sAs', 'sBs', 'sCs', 'alph', 'beta', 'gamm']
	}, 
	'vbaimft': {
		'n': 9,
		'p': ['tau', 'sAs', 'sBs', 'sCs', 'Ea', 'alph', 'beta', 'gamm']
	}, 
	
	'bgaf': {
		'n': 4,
		'p': ['tau', 'sAs', 'sBs', 'sGs']
	}, 
	'bgaft': {
		'n': 6,
		'p': ['tau', 'sAs', 'sBs', 'sGs', 'Ea']
	}, 
	'vbgaf': {
		'n': 7,
		'p': ['tau', 'sAs', 'sBs', 'sGs', 'alph', 'beta', 'gamm']
	}, 
	'vbgaft': {
		'n': 9,
		'p': ['tau', 'sAs', 'sBs', 'sGs', 'Ea', 'alph', 'beta', 'gamm']
	}, 
	
	
	'uvegaf': {
		'n': 10,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sGs', 'S2f', 'S2uf', 'alph', 'beta', 'gamm']
	}, 
	'uvegaft': {
		'n': 12,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sGs', 'S2f', 'Eas', 'Eaf', 'S2uf', 'alph', 'beta', 'gamm']
	}, 
	'uvgaf': {
		'n': 12,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sGs', 'sAf', 'sBf', 'sGf', 'S2uf', 'alph', 'beta', 'gamm']
	}, 
	'uvgaft': {
		'n': 14,
		'p': ['taus', 'tauf', 'sAs', 'sBs', 'sGs', 'sAf', 'sBf', 'sGf', 'Eas', 'Eaf', 'S2uf', 'alph', 'beta', 'gamm']
	}, 
}
