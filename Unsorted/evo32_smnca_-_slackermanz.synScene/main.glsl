



			//******** BuffA Code Begins ********



//	----    ----    ----    ----    ----    ----    ----    ----

//	NAME:Evo32 Selective MNCA

//	TYPE:Selective Multiple Neighbourhood Cellular Automata

//	RULE:Pattern coordinates hardcoded into 'demo' array

//	----    ----    ----    ----    ----    ----    ----    ----

//	Specification: SlackCA_Specification_2021_02_23:

//		https://mega.nz/file/yRliVDJT#6CUlcGma4DpfXI4S8j0VoUi8Vju0vwRXwVI4klyiNXg

//	Adapted for Shadertoy. FPS compared to the C++/Vulkan application is about 20% of maximum

//

//	Text file containing other pattern coordinates:

//		https://mega.nz/file/uY0GjSjJ#-TvjklejZBh3O5DfqtLXFtjZaoVetqHPgotdYLH5xoQ

//	Visualisation of a small subset of rules recorded by Softology:

//		https://www.youtube.com/watch?v=LtdKNso0DwE

//

//	----    ----    ----    ----    ----    ----    ----    ----

//  Shader developed by Slackermanz

//

//  Info/Code:

//  ﻿ - Website: https://slackermanz.com

//  ﻿ - Github: https://github.com/Slackermanz

//  ﻿ - Shadertoy: https://www.shadertoy.com/user/SlackermanzCA

//  ﻿ - Discord: https://discord.gg/hqRzg74kKT

//  

//  Socials:

//  ﻿ - Discord DM: Slackermanz#3405

//  ﻿ - Reddit DM: https://old.reddit.com/user/slackermanz

//  ﻿ - Twitter: https://twitter.com/slackermanz

//  ﻿ - YouTube: https://www.youtube.com/c/slackermanz

//  ﻿ - Older YT: https://www.youtube.com/channel/UCZD4RoffXIDoEARW5aGkEbg

//  ﻿ - TikTok: https://www.tiktok.com/@slackermanz

//  

//  Communities:

//  ﻿ - Reddit: https://old.reddit.com/r/cellular_automata

//  ﻿ - Artificial Life: https://discord.gg/7qvBBVca7u

//  ﻿ - Emergence: https://discord.com/invite/J3phjtD

//  ﻿ - ConwayLifeLounge: https://discord.gg/BCuYCEn

//	----    ----    ----    ----    ----    ----    ----    ----



#define tex0_in (BuffA)

#define tex_size (RENDERSIZE.xy)



precision highp int;

precision highp float;



float gdv(float x, float y, int v, float div) {

//	Get Div Value: Return the value of a specified pixel

//		x, y : 	Relative integer-spaced coordinates to origin [ 0.0, 0.0 ]

//		v	 :	Colour channel [ 0, 1, 2 ]

//		div	 :	Integer-spaced number of toroidal divisions of the surface/medium

	vec4	fc 		= gl_FragCoord;

	vec2 	dm 		= vec2(RENDERSIZE.x, RENDERSIZE.y);

    float 	divx	= dm[0] / div;

	float 	divy	= dm[1] / div;

	float 	pxo 	= 1.0 / dm[0];

	float 	pyo 	= 1.0 / dm[1];

	float 	fcxo 	= fc[0] + x;

	float 	fcyo 	= fc[1] + y;

	float 	fcx 	= (mod(fcxo,divx) + floor(fc[0]/divx)*divx ) * pxo;

	float 	fcy 	= (mod(fcyo,divy) + floor(fc[1]/divy)*divy ) * pyo;

    vec4 	pxdata 	= texture( tex0_in, vec2(fcx, fcy) );

	return pxdata[v];

}



vec3 nhd( vec2 nbhd, vec2 ofst, float psn, float thr, int col, float div ) {

//	Neighbourhood: Return information about the specified group of pixels

	float dist 		= 0.0;

	float cval 		= 0.0;

	float c_total 	= 0.0;

	float c_valid 	= 0.0;

	float c_value 	= 0.0;

	for(float i = -nbhd[0]; i <= nbhd[0]; i += 1.0) {

		for(float j = -nbhd[0]; j <= nbhd[0]; j += 1.0) {

			dist = round(sqrt(i*i+j*j));

			if( dist <= nbhd[0] && dist > nbhd[1] && dist != 0.0 ) {

				cval = gdv(i+ofst[0],j+ofst[1],col,div);

				c_total += psn;

				if( cval > thr ) {

					c_valid += psn;

					cval = psn * cval;

					c_value += cval-fract(cval); } } } } 

	return vec3( c_value, c_valid, c_total );

}



float get_xc(float x, float y, float xmod) {

//	Used to reseed the surface with noise

	float sq = sqrt(mod(x*y+y, xmod)) / sqrt(xmod);

	float xc = mod((x*x)+(y*y), xmod) / xmod;

	return clamp((sq+xc)*0.5, 0.0, 1.0);

}

float shuffle(float x, float y, float xmod, float val) {

//	Used to reseed the surface with noise

	val = val * mod( x*y + x, xmod );

	return (val-floor(val));

}

float get_xcn(float x, float y, float xm0, float xm1, float ox, float oy) {

//	Used to reseed the surface with noise

	float  xc = get_xc(x+ox, y+oy, xm0);

	return shuffle(x+ox, y+oy, xm1, xc);

}

float get_lump(float x, float y, float nhsz, float xm0, float xm1) {

//	Used to reseed the surface with noise

	float 	nhsz_c 	= 0.0;

	float 	xcn 	= 0.0;

	float 	nh_val 	= 0.0;

	for(float i = -nhsz; i <= nhsz; i += 1.0) {

		for(float j = -nhsz; j <= nhsz; j += 1.0) {

			nh_val = round(sqrt(i*i+j*j));

			if(nh_val <= nhsz) {

				xcn = xcn + get_xcn(x, y, xm0, xm1, i, j);

				nhsz_c = nhsz_c + 1.0; } } }

	float 	xcnf 	= ( xcn / nhsz_c );

	float 	xcaf	= xcnf;

	for(float i = 0.0; i <= nhsz; i += 1.0) {

			xcaf 	= clamp((xcnf*xcaf + xcnf*xcaf) * (xcnf+xcnf), 0.0, 1.0); }

	return xcaf;

}

float reseed() {

//	Used to reseed the surface with noise

	vec4	fc = gl_FragCoord;

	float 	r0 = get_lump(fc[0], fc[1], 2.0, 19.0 + mod(vec4(2019.0, 1.0, 1.0, TIME)[3],17.0), 23.0 + mod(vec4(2019.0, 1.0, 1.0, TIME)[3],43.0));

	float 	r1 = get_lump(fc[0], fc[1], 24.0, 13.0 + mod(vec4(2019.0, 1.0, 1.0, TIME)[3],29.0), 17.0 + mod(vec4(2019.0, 1.0, 1.0, TIME)[3],31.0));

	float 	r2 = get_lump(fc[0], fc[1], 8.0, 13.0 + mod(vec4(2019.0, 1.0, 1.0, TIME)[3],11.0), 51.0 + mod(vec4(2019.0, 1.0, 1.0, TIME)[3],37.0));

	return clamp((r0+r1)-r2,0.0,1.0);

}



vec4 renderPassA() {

	vec4 fragColor = vec4(0.0);

	vec2 fragCoord = _xy;



//	----    ----    ----    ----    ----    ----    ----    ----

//	Initilisation

//	----    ----    ----    ----    ----    ----    ----    ----

			vec4	fc 		= gl_FragCoord;							// 	Frag Coords

			vec2 	dm 		= vec2(RENDERSIZE.x, RENDERSIZE.y);	//	Surface Dimensions

	const 	int 	VMX 	= 32;									//	Evolution Variables

			float 	psn		= 250.0;								//	Precision

			float 	mnp		= 0.004;								//	Minimum Precision Value : (1.0 / psn);

			float	div		= 2.0;									//	Toroidal Surface Divisions

			float 	divi 	= floor((fc[0]*div)/(dm[0]))		

							+ floor((fc[1]*div)/(dm[1]))*div;		//	Division Index

	float 	dspace  = (divi+1.0)/(div*div);

			dspace  = (div == 1.0) ? 1.0 : dspace;					//	Division Weight

	vec3 	col 	= vec3( 0.0, 0.0, 0.0 );						//	Final colour value output container



//	Pattern coordinates hardcoded for demonstration purposes

//	Text file containing other pattern coordinates:

//		https://mega.nz/file/uY0GjSjJ#-TvjklejZBh3O5DfqtLXFtjZaoVetqHPgotdYLH5xoQ

	float[VMX] 	demo_0 	= float[VMX]

	(	0.264000,	0.293000,	0.156000,	0.287000,	

		0.411000,	0.305000,	0.200000,	0.093000,	

		0.467000,	0.345000,	0.127000,	-0.016000,	

		-0.108000,	0.079000,	0.098000,	0.370000,	

		0.209000,	0.030000,	0.222000,	0.144000,	

		0.187000,	0.199000,	0.294000,	0.194000,	

		0.277000,	0.170000,	0.178000,	0.007000,	

		0.261000,	0.335000,	-0.093000,	0.138000	);

        

	float[VMX] 	demo_1 	= float[VMX]

	(	0.352000,	0.236000,	0.313000,	0.197000,	

		0.140000,	0.443000,	0.252000,	0.394000,	

		0.450000,	0.567000,	0.233000,	0.387000,	

		0.273000,	0.489000,	0.190000,	0.374000,	

		0.215000,	0.241000,	0.264000,	0.387000,	

		0.308000,	0.373000,	0.263000,	0.499000,	

		0.225000,	0.252000,	0.168000,	0.344000,	

		0.179000,	0.324000,	0.391000,	0.327000	);

        

	float[VMX] 	demo_2 	= float[VMX]

	(	0.534000,	-0.129000,	0.041000,	0.683000,	

		0.228000,	0.562000,	-0.097000,	0.340000,	

		-0.043000,	0.209000,	-0.136000,	0.746000,	

		0.663000,	0.222000,	0.768000,	0.299000,	

		0.354000,	0.449000,	0.172000,	0.263000,	

		0.074000,	0.760000,	-0.020000,	0.060000,	

		0.133000,	0.198000,	0.200000,	0.516000,	

		0.320000,	0.429000,	0.378000,	0.710000	);

        

	float[VMX] 	demo_3 	= float[VMX]

	(	0.378000,	0.228000,	0.196000,	-0.046000,	

		0.156000,	-0.055000,	0.017000,	0.055000,	

		0.065000,	0.255000,	0.115000,	0.395000,	

		-0.041000,	0.029000,	0.000000,	0.125000,	

		0.352000,	0.346000,	0.178000,	0.305000,	

		0.020000,	0.237000,	0.397000,	0.311000,	

		0.215000,	0.048000,	0.344000,	-0.009000,	

		0.359000,	0.244000,	0.022000,	0.174000	);



//	Division Weighted V number

	float[VMX] 	dvmd;

    if(divi == 0.0) {

		for(int i = 0; i < VMX; i++) { dvmd[i] = demo_0[i]; } }

    if(divi == 1.0) {

		for(int i = 0; i < VMX; i++) { dvmd[i] = demo_1[i]; } }

    if(divi == 2.0) {

		for(int i = 0; i < VMX; i++) { dvmd[i] = demo_2[i]; } }

    if(divi == 3.0) {

		for(int i = 0; i < VMX; i++) { dvmd[i] = demo_3[i]; } }

    

//	----    ----    ----    ----    ----    ----    ----    ----

//	RULE:Evo32 Selective MNCA

//	TYPE:Selective Multiple Neighbourhood Cellular Automata

//	----    ----    ----    ----    ----    ----    ----    ----

//	Rule Initilisation



//	Get the reference frame's origin pixel values

	float	res_r	= gdv( 0.0, 0.0, 0, div );

	float	res_g	= gdv( 0.0, 0.0, 1, div );

	float	res_b	= gdv( 0.0, 0.0, 2, div );



//	Intended rate of change

	float 	s 		= mnp * 12.0;



//	----    ----    ----    ----    ----    ----    ----    ----

//	STAGE:MNCA

//		DOMAIN: 	Totalistic Multiple Neighbourhood Continuous

//		REQUIRES:	Conditional Range

//		UPDATE:		Additive

//		VALUE:		Origin

//		BLUR:		Relative

//		RESULT:		Multiple



//	Number of Individual Neighbourhoods

	const 	int 		nhds	= 8;

//	Container for Neighbourhood Totalistic values

			float[nhds] nhdt;



//	Number of MNCA Groups

	const 	int 		sets	= 4;

//	Container for STAGE:MNCA results

			float[sets] rslt;

						rslt[0] = res_r;

						rslt[1] = res_r;

						rslt[2] = res_r;

						rslt[3] = res_r;



//	Define and assess the Individual Neighbourhoods

	vec3 	n00r 	= nhd( vec2( 1.0, 	0.0 ), vec2( 0.0, 0.0 ), psn, 0.0, 0, div );

	vec3 	n01r 	= nhd( vec2( 3.0, 	0.0 ), vec2( 0.0, 0.0 ), psn, 0.0, 0, div );

	vec3 	n10r 	= nhd( vec2( 2.0, 	0.0 ), vec2( 0.0, 0.0 ), psn, 0.0, 0, div );

	vec3 	n11r 	= nhd( vec2( 5.0, 	3.0 ), vec2( 0.0, 0.0 ), psn, 0.0, 0, div );

	vec3 	n20r 	= nhd( vec2( 5.0, 	2.0 ), vec2( 0.0, 0.0 ), psn, 0.0, 0, div );

	vec3 	n21r 	= nhd( vec2( 9.0, 	7.0 ), vec2( 0.0, 0.0 ), psn, 0.0, 0, div );

	vec3 	n30r 	= nhd( vec2( 4.0, 	2.0 ), vec2( 0.0, 0.0 ), psn, 0.0, 0, div );

	vec3 	n31r 	= nhd( vec2( 12.0, 	9.0 ), vec2( 0.0, 0.0 ), psn, 0.0, 0, div );



//	Get the Totalistic value of each Individual Neighbourhood

	nhdt[0] = n00r[0] / n00r[2];

	nhdt[1] = n01r[0] / n01r[2];

	nhdt[2] = n10r[0] / n10r[2];

	nhdt[3] = n11r[0] / n11r[2];

	nhdt[4] = n20r[0] / n20r[2];

	nhdt[5] = n21r[0] / n21r[2];

	nhdt[6] = n30r[0] / n30r[2];

	nhdt[7] = n31r[0] / n31r[2];

    



//	UPDATE:Additive the VALUE:Origin according to REQUIRES:Conditional Range

	if(	nhdt[0] >= dvmd[0] 	&& nhdt[0] <= dvmd[1] 	) { rslt[0] += s; }

	if(	nhdt[0] >= dvmd[2] 	&& nhdt[0] <= dvmd[3] 	) { rslt[0] -= s; }

	if( nhdt[1] >= dvmd[4] 	&& nhdt[1] <= dvmd[5] 	) { rslt[0] += s; }

	if( nhdt[1] >= dvmd[6] 	&& nhdt[1] <= dvmd[7] 	) { rslt[0] -= s; }



	if( nhdt[2] >= dvmd[8] 	&& nhdt[2] <= dvmd[9] 	) { rslt[1] += s; }

	if( nhdt[2] >= dvmd[10] && nhdt[2] <= dvmd[11] 	) { rslt[1] -= s; }

	if( nhdt[3] >= dvmd[12] && nhdt[3] <= dvmd[13] 	) { rslt[1] += s; }

	if( nhdt[3] >= dvmd[14] && nhdt[3] <= dvmd[15] 	) { rslt[1] -= s; }



	if(	nhdt[4] >= dvmd[16] && nhdt[4] <= dvmd[17] 	) { rslt[2] += s; }

	if(	nhdt[4] >= dvmd[18] && nhdt[4] <= dvmd[19] 	) { rslt[2] -= s; }

	if( nhdt[5] >= dvmd[20] && nhdt[5] <= dvmd[21] 	) { rslt[2] += s; }

	if( nhdt[5] >= dvmd[22] && nhdt[5] <= dvmd[23] 	) { rslt[2] -= s; }



	if( nhdt[6] >= dvmd[24] && nhdt[6] <= dvmd[25] 	) { rslt[3] += s; }

	if( nhdt[6] >= dvmd[26] && nhdt[6] <= dvmd[27] 	) { rslt[3] -= s; }

	if( nhdt[7] >= dvmd[28] && nhdt[7] <= dvmd[29] 	) { rslt[3] += s; }

	if( nhdt[7] >= dvmd[30] && nhdt[7] <= dvmd[31] 	) { rslt[3] -= s; }



//	Apply a BLUR:Relative to the result

	rslt[0] = (rslt[0] + nhdt[0] * s * 0.5 + nhdt[1] * s * 0.5) / (1.0 + s * 1.0 );

	rslt[1] = (rslt[1] + nhdt[2] * s * 0.5 + nhdt[3] * s * 0.5) / (1.0 + s * 1.0 );

	rslt[2] = (rslt[2] + nhdt[4] * s * 0.5 + nhdt[5] * s * 0.5) / (1.0 + s * 1.0 );

	rslt[3] = (rslt[3] + nhdt[6] * s * 0.5 + nhdt[7] * s * 0.5) / (1.0 + s * 1.0 );



//	----    ----    ----    ----    ----    ----    ----    ----

//	STAGE:Variance

//		DOMAIN: 	MNCA

//		REQUIRES:	Unconditional

//		UPDATE:		Subtract

//		VALUE:		Origin

//		BLUR:		Specific

//		RESULT:		Multiple



//	Container for STAGE:Variance results

	float[sets] variance;



//	UPDATE:Subtract the REQUIRES:Previous value

	for(int i = 0; i < sets; i++) { 

		{ variance[i] = res_r - rslt[i]; } }



//	----    ----    ----    ----    ----    ----    ----    ----

//	STAGE:Output

//		DOMAIN: 	MNCA, Variance

//		REQUIRES:	MaximumABS_Match

//		UPDATE:		Select

//		VALUE:		Domain[0]

//		BLUR:		Specific

//		RESULT:		Single



//	Index of an element in DOMAIN:Variance

	int von = 0;



//	Get the index of the element in DOMAIN:Variance that meets REQUIRES:MaximumABS_Match

	for( int i = 0; i < sets; i++ ) { if( abs(variance[von]) < abs(variance[i]) ) { von = i; } }



//	UPDATE:Select the DOMAIN:MNCA value with the REQUIRES:MaximumABS_Match index

	float maxvar = rslt[von];



//	Output that value

	res_r = maxvar;



//	----    ----    ----    ----    ----    ----    ----    ----

//	Presentation Filtering

//	----    ----    ----    ----    ----    ----    ----    ----

	vec3 	n0g 	= nhd( vec2( 1.0, 0.0 ), vec2( 0.0, 0.0 ), psn, 0.0, 1, div );

	vec3 	n0b 	= nhd( vec2( 2.0, 0.0 ), vec2( 0.0, 0.0 ), psn, 0.0, 2, div );

	float 	n0gw 	= n0g[0] / n0g[2];

	float 	n0bw 	= n0b[0] / n0b[2];

	res_g = ( res_g + n0gw * s * 2.0 + res_r * s * 2.0 ) / (1.0 + s * 4.0);

	res_b = ( res_b + n0bw * s * 1.0 + res_r * s * 1.0 ) / (1.0 + s * 2.0);



//	----    ----    ----    ----    ----    ----    ----    ----

//	Fragment Shader Output

//	----    ----    ----    ----    ----    ----    ----    ----

    

    if (_mouse.z > 0. && length(_mouse.xy - fragCoord) < 32.0) { res_r = mod(float(FRAMECOUNT),2.0); }

    if (FRAMECOUNT == 0) { res_r = reseed(); }

    fragColor=vec4(clamp(res_r,0.0,1.0),clamp(res_g,0.0,1.0),clamp(res_b,0.0,1.0),1.0);

	return fragColor; 

 } 





//	----    ----    ----    ----    ----    ----    ----    ----

//	NAME:Evo32 Selective MNCA

//	TYPE:Selective Multiple Neighbourhood Cellular Automata

//	RULE:Pattern coordinates hardcoded into 'demo' array

//	----    ----    ----    ----    ----    ----    ----    ----

//	Specification: SlackCA_Specification_2021_02_23:

//		https://mega.nz/file/yRliVDJT#6CUlcGma4DpfXI4S8j0VoUi8Vju0vwRXwVI4klyiNXg

//	Adapted for Shadertoy. FPS compared to the C++/Vulkan application is about 20% of maximum

//

//	Text file containing other pattern coordinates:

//		https://mega.nz/file/uY0GjSjJ#-TvjklejZBh3O5DfqtLXFtjZaoVetqHPgotdYLH5xoQ

//	Visualisation of a small subset of rules recorded by Softology:

//		https://www.youtube.com/watch?v=LtdKNso0DwE

//

//	----    ----    ----    ----    ----    ----    ----    ----

//  Shader developed by Slackermanz

//

//  Info/Code:

//  ﻿ - Website: https://slackermanz.com

//  ﻿ - Github: https://github.com/Slackermanz

//  ﻿ - Shadertoy: https://www.shadertoy.com/user/SlackermanzCA

//  ﻿ - Discord: https://discord.gg/hqRzg74kKT

//  

//  Socials:

//  ﻿ - Discord DM: Slackermanz#3405

//  ﻿ - Reddit DM: https://old.reddit.com/user/slackermanz

//  ﻿ - Twitter: https://twitter.com/slackermanz

//  ﻿ - YouTube: https://www.youtube.com/c/slackermanz

//  ﻿ - Older YT: https://www.youtube.com/channel/UCZD4RoffXIDoEARW5aGkEbg

//  ﻿ - TikTok: https://www.tiktok.com/@slackermanz

//  

//  Communities:

//  ﻿ - Reddit: https://old.reddit.com/r/cellular_automata

//  ﻿ - Artificial Life: https://discord.gg/7qvBBVca7u

//  ﻿ - Emergence: https://discord.com/invite/J3phjtD

//  ﻿ - ConwayLifeLounge: https://discord.gg/BCuYCEn

//	----    ----    ----    ----    ----    ----    ----    ----



precision highp int;

precision highp float;



vec4 renderMainImage() {

	vec4 fragColor = vec4(0.0);

	vec2 fragCoord = _xy;



    // Normalized pixel coordinates (from 0 to 1)

    vec2 uv = fragCoord/RENDERSIZE.xy;



    vec4 bufA = texture(BuffA, uv);



    // Output to screen

    fragColor = bufA;

	return fragColor; 

 } 





vec4 renderMain(){

	if(PASSINDEX == 0){

		return renderPassA();

	}

	if(PASSINDEX == 1){

		return renderMainImage();

	}

}