



			//******** BuffA Code Begins ********


vec4 mediaEdges = texture(media_pass_fx, _uv)*media_impact*0.125;
vec4 mixedEdgeCol = mix(_loadMediaAsMask(), mediaEdges, .25+0.75*edge_mix);
// vec4 mediaMixer = mix(texture(media_pass, _uv), texture(media_pass_fx, _uv), 1.0-media_color_mix);
float mediahits_mix =  mix(media_impact, media_impact* syn_MidLevel*0.1*syn_Intensity +syn_Presence*syn_Level*0.125, media_hits);
bool mediaOn = _exists(syn_UserImage);
float media_lum = PI*sin(length(mixedEdgeCol*mediahits_mix));
float m_lum = mix(0, media_lum, int(mediaOn));
// float paintsize = paint_size * 40.;


#define D(COORD) texture(BuffD,(COORD)/RENDERSIZE.xy)
float growthFactor =0.25+2*normalize(pow((syn_BassLevel*0.5)+(syn_MidLevel*0.25)+(syn_Intensity*0.25), 2.0));

vec4 image = vec4(0.);


vec4 renderPassA() {

	vec4 Q = vec4(0.0);

	vec2 U = _xy;
    vec3 edges = _edgeDetectSobel(syn_UserImage).rgb;
    U += vec2(_noise(cos(TIME)), _noise(sin(TIME)));
    U += _fbm(_uvc*PI)*Warp;
    U -= _uvc*PI*Zoom*(1.0+pow(syn_MidLevel*0.5+syn_BassLevel*0.5, 2.));
    vec3 logoCol = _loadMedia().rgb;

     
    image.rgb = mix(edges.rgb, logoCol.rgb, .75-syn_Level);




    U += MoveXY*(1.0+pow(syn_MidLevel*0.5+syn_BassLevel*0.5, 2.));

    U += Stretch*_uvc*2.;

    vec2 X = 0.5*RENDERSIZE.xy;

    if (_mouse.z>0.) X = _mouse.xy;

    U-=X;

    float a = .001*sin(.01*smoothTimeC)/(1.+length(U-0.5*RENDERSIZE.xy)/RENDERSIZE.y);

    
    U *= mat2(cos(a),-sin(a),sin(a),cos(a));

    U *= .999;

    U+=X;

    //U += (_loadUserImageAsMask().rb*0.1*syn_Intensity*_uvc);



    Q  =  D(U);

    // Neighborhood :

    vec4 pX  =  D(U + vec2(1,0));
    vec4 pY  =  D(U + vec2(0,1));
    vec4 nX  =  D(U - vec2(1,0));
    vec4 nY  =  D(U - vec2(0,1));

    vec4 m = 0.25*(pX+nX+pY+nY);

        m = mix(image, m, 0.975-syn_BassLevel*0.1);

    //float b = .3;

        float b = .3+Size;



    Q += (1.-b)*(0.25*vec4(

       

        (pX.z-nX.z) - 5.*Q.z*(pY.z-nY.z),

       	(pY.z-nY.z) + 5.*Q.z*(pX.z-nX.z),

        -pX.x+nX.x-pY.y+nY.y,1

    )- Q);


    image = mix(image, Q, 0.125);

    Q = mix(Q,m,b);


    
    Q = mix(Q, image, 0.0125);

    if (length(Q.xy)>0.) Q.xy = normalize(Q.xy);

    

    if(FRAMECOUNT <= 1 || Reset != 0.) Q = sin(.1*U.xxxx +_uvc.xxxx);

     

	return Q; 

 } 





			//******** BuffB Code Begins ********



#define A(COORD) texture(BuffA,(COORD)/RENDERSIZE.xy)



vec4 renderPassB() {

	vec4 Q = vec4(0.0);

	vec2 U = _xy;

    

    U += Stretch*_uvc*2.;

    U += MoveXY*(1.0+pow(syn_MidLevel*0.5+syn_BassLevel*0.5, 2.));



    vec2 X = 0.5*RENDERSIZE.xy;

    if (_mouse.z>0.) X = _mouse.xy;

    U-=X;

    float a = .001*sin(.01*smoothTimeC)/(1.+length(U-0.5*RENDERSIZE.xy)/RENDERSIZE.y);

    U *= mat2(cos(a),-sin(a),sin(a),cos(a));

    U *= .999;

    U+=X;

    Q  =  A(U);

    // Neighborhood :

    vec4 pX  =  A(U + vec2(1,0));
    vec4 pY  =  A(U + vec2(0,1));
    vec4 nX  =  A(U - vec2(1,0));
    vec4 nY  =  A(U - vec2(0,1));

    vec4 m = 0.25*(pX+nX+pY+nY);
    float b = .29;

    //b -= syn_BassLevel*0.075;

    Q += (1.-b)*(0.25*vec4(

       

        (pX.z-nX.z) - 5.*Q.z*(pY.z-nY.z),

       	(pY.z-nY.z) + 5.*Q.z*(pX.z-nX.z),

        -pX.x+nX.x-pY.y+nY.y,1

    )- Q);




    
    Q = mix(Q*(0.75+0.3*syn_MidLevel+0.3*syn_BassLevel),m,b);

  

    if (length(Q.xy)>0.) Q.xy = normalize(Q.xy);

    

    if(FRAMECOUNT <= 1 || Reset != 0.) Q = sin(.1*U.xxxx);

    

     

	return Q; 

 } 





			//******** BuffC Code Begins ********



#define B(COORD) texture(BuffB,(COORD)/RENDERSIZE.xy)



vec4 renderPassC() {

	vec4 Q = vec4(0.0);

	vec2 U = _xy;

    U += MoveXY*(1.0+pow(syn_MidLevel*0.5+syn_BassLevel*0.5, 2.));

    U -= _uvc*PI*Zoom*(1.0+pow(syn_MidLevel*0.5+syn_BassLevel*0.5, 2.));



    vec2 X = 0.5*RENDERSIZE.xy;

    if (_mouse.z>0.) X = _mouse.xy;

    U-=X;

    float a = .001*sin(.01*smoothTimeC)/(1.+length(U-0.5*RENDERSIZE.xy)/RENDERSIZE.y);

    U *= mat2(cos(a),-sin(a),sin(a),cos(a));

    U *= .999;

    U+=X;

    Q  =  B(U);

    // Neighborhood :

    vec4 pX  =  B(U + vec2(1,0));

    vec4 pY  =  B(U + vec2(0,1));

    vec4 nX  =  B(U - vec2(1,0));

    vec4 nY  =  B(U - vec2(0,1));

    vec4 m = 0.25*(pX+nX+pY+nY);

    float b = .3;

    Q += (1.-b)*(0.125*(1.0+syn_MidLevel)*vec4(

       

        (pX.z-nX.z) - 5.*Q.z*(pY.z-nY.z),

       	(pY.z-nY.z) + 5.*Q.z*(pX.z-nX.z),

        -pX.x+nX.x-pY.y+nY.y,1

    )- Q);



    

    Q = mix(Q,m,b);

    

    if (length(Q.xy)>0.) Q.xy = normalize(Q.xy);

    

    if(FRAMECOUNT <= 1 || Reset != 0.) Q = sin(.1*U.xxxx);

    

     

	return Q; 

 } 





			//******** BuffD Code Begins ********



#define C(COORD) texture(BuffC,(COORD)/RENDERSIZE.xy)



vec4 renderPassD() {

	vec4 Q = vec4(0.0);

	vec2 U = _xy+2.*_uvc/PI;



    vec2 X = 0.5*RENDERSIZE.xy;

    if (_mouse.z>0.) X = _mouse.xy;

    U-=X;

    float a = .001*sin(.01*smoothTimeC)/(1.+length(U-0.5*RENDERSIZE.xy)/RENDERSIZE.y);

    U *= mat2(cos(a),-sin(a),sin(a),cos(a));

    U *= .999;

    U+=X;

    

    Q  =  C(U);

    

    // Neighborhood :

    vec4 pX  =  C(U + vec2(1,0)+media_lum);

    vec4 pY  =  C(U + vec2(0,1)+media_lum);

    vec4 nX  =  C(U - vec2(1,0)+media_lum);

    vec4 nY  =  C(U - vec2(0,1)+media_lum);

    vec4 m = 0.25*(pX+nX+pY+nY);

    float b = .3;

    Q += (1.-b)*(0.25*vec4(

       

        (pX.z-nX.z) - 5.*Q.z*(pY.z-nY.z),

       	(pY.z-nY.z) + 5.*Q.z*(pX.z-nX.z),

        -pX.x+nX.x-pY.y+nY.y,1

    )- Q);



    

    Q = mix(Q,m,b);

    

    if (length(Q.xy)>0.) Q.xy = normalize(Q.xy);

    

    if(FRAMECOUNT <= 1 || Reset != 0.) Q = sin(.1*U.xxxx);

    

     

	return Q; 

 } 





#define A(COORD) texture(BuffA,(COORD)/RENDERSIZE.xy)

float ln (vec3 p, vec3 a, vec3 b) {return length(p-a-(b-a)*min(dot(p-a,b-a),0.)/dot(b-a,b-a));}

vec4 renderMainImage() {

	vec4 Q = vec4(0.0);

	vec2 U = _xy;



    vec2 X = vec2(0.5)+2.*_uvc/PI;

    if (_mouse.z>0.) X = _mouse.xy/RENDERSIZE.xy;

    Q  =  A(U);

    vec4 pX  =  A(U + vec2(1,0));

    vec4 pY  =  A(U + vec2(0,1));

    vec4 nX  =  A(U - vec2(1,0));

    vec4 nY  =  A(U - vec2(0,1));

    vec3 n = normalize(vec3(pX.z-nX.z,pY.z-nY.z,2));

    vec3 r = reflect(n,vec3(0,0,-1));

    Q = 0.125+0.5*sin((2.5)*vec4(1,2,3,4)*(Q.z));

    float d = ln(vec3(X+.5,8)*RENDERSIZE.xyy,

                 vec3(U,0),vec3(U,0)+r)/RENDERSIZE.y;

    Q *= exp(-(.5)*d*d)+0.5*pow(syn_HighLevel*0.5+0.5*syn_Intensity, 2.);

	return Q; 

 } 


vec4 mediaPass() {
    vec4 media = vec4(0.);
    vec2 U = _xy;
    media = _loadMedia();
    return media;
}
vec4 mediaPassFX() {
    vec4 media_fx = _edgeDetectSobel(media_pass);
    return media_fx;
} 




vec4 renderMain(){

	if(PASSINDEX == 0){

		return renderPassA();

	}

	if(PASSINDEX == 1){

		return renderPassB();

	}

	if(PASSINDEX == 2){

		return renderPassC();

	}

	if(PASSINDEX == 3){

		return renderPassD();

	}

	if(PASSINDEX == 4){

		return mediaPass();

	}
	if(PASSINDEX == 5){

		return mediaPassFX();

	}
	if(PASSINDEX == 6){

		return renderMainImage();

	}

}