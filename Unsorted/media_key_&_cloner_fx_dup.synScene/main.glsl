#define TAU             (2.0*PI)
#define ROT(a)          mat2(cos(a), sin(a), -sin(a), cos(a))
#define NSIN(x)         (0.5+0.5*sin(x))
#define LESS(a,b,c)     mix(a,b,step(0.,c))
#define SABS(x,k)       LESS((.5/(k))*(x)*(x)+(k)*.5,abs(x),abs(x)-(k))

//helper functions
float loop(float x){return abs(fract(x*0.5+0.5)*2.0-1.0);}

vec2 RGBtoYUV(vec3 rgb) {
  return vec2(
    rgb.r * -0.169 + rgb.g * -0.331 + rgb.b *  0.5    + 0.5,
    rgb.r *  0.5   + rgb.g * -0.419 + rgb.b * -0.081  + 0.5
  );
}

//from mrange: https://www.shadertoy.com/view/WlcfRS
vec2 toPolar(vec2 p) {
  return vec2(length(p), atan(p.y, p.x));
}

vec2 toRect(vec2 p) {
  return vec2(p.x*cos(p.y), p.x*sin(p.y));
}

float modMirror1(inout float p, float size) {
  float halfsize = size*0.5;
  float c = floor((p + halfsize)/size);
  p = mod(p + halfsize,size) - halfsize;
  p *= mod(c, 2.0)*2.0 - 1.0;
  return c;
}

vec2 smoothKaleidoscope(vec2 p, float sm, float rep) {
  vec2 hp = p;
  vec2 hpp = toPolar(hp);
  float rn = modMirror1(hpp.y, TAU/rep);
  float sa = PI/rep - SABS(PI/rep - abs(hpp.y), sm);
  hpp.y = sign(hpp.y)*(sa);
  hp = toRect(hpp);
  p = hp;
//   return rn;
  return hp;
}

vec4 getContentMask(vec4 fragColor) {
    vec4 contentMask = fragColor;
    
    /////////////////////
    //color
    ///////////////////
    float tolerance = 0.46-matte_tolerance*0.45;
    vec3 matteColor = mix(color_key, texture(keyBuffer, vec2(0.0)).rgb, clamp(keying_style - 2., 0.0, 1.0));
    float value = distance(RGBtoYUV(fragColor.rgb), RGBtoYUV(matteColor));
    float valueFeather = distance(RGBtoYUV(fragColor.rgb), RGBtoYUV(matteColor + feather_edge));
    float colorMask = value - tolerance;
    float colorKey = pow(clamp(colorMask / 0.08, 0., 1.), matte_choke*2.);
    float offsetColorKey = pow(clamp(colorMask / 0.17, 0., 1.), matte_choke*2.);
    float offsetColorKey2 = pow(clamp(colorMask / 0.23, 0., 1.), matte_choke*2.);
    colorKey = colorKey - (colorKey - offsetColorKey)*feather_edge;
    
            
    // For spill
    float desat = clamp(contentMask.r * 0.01 + contentMask.g * 0.01 + contentMask.b * 0.01, 0., 1.);
    contentMask.rgb = mix(vec3(desat, desat, desat)*colorKey, contentMask.rgb, offsetColorKey2);
    // contentMask.rgb = mix(vec3(desat, desat, desat), contentMask.rgb, offsetColorKey*step(1.5, keying_style)*step(keying_style, 3.5));
        
        
    ///////////////
    //luma
    //////////////
    float luma = (length(fragColor.rgb - (1.0 - matte_tolerance)));
    
    float lumaKey = smoothstep(-(1.0-matte_choke), feather_edge+0.01,luma-tolerance);
    float lumaGate = smoothstep(0.0, 0.6, step(matte_choke, luma) + smoothstep(matte_choke*0.7, matte_choke, luma + feather_edge*0.2)*0.9*(dot(fragColor.rgb, vec3(1.))/3.));
    lumaKey = mix(lumaKey, lumaGate, clamp(keying_style, 0.0, 1.0));

    
    ////////////////
    // optical flow
    //////////////////////
    vec3 opticalFlow = texture(opticalFlowBuffer, _uv).rgb;
    float opticalFlowLuma = dot(opticalFlow.rgb, vec3(1.))/3.*matte_choke;
    float opticalFlowKey = smoothstep(0.0,0.3 + feather_edge+0.01,opticalFlowLuma);
    
    contentMask.a = mix(lumaKey, mix(colorKey, opticalFlowKey, clamp(keying_style - 3., 0.0, 1.0)), clamp(keying_style - 1.0, 0.0, 1.0));

    return contentMask;
}

float getVignette() {
    return smoothstep(0.0, vignette_softness+0.01, (1.0 - vignette_amount) - (1.0-vignette_softness)*0.25 + (1.0 - length(mix(_uvc*0.97, vec2(abs(2.0*_uv-1.0)), vignette_roundness))));
}


vec4 renderBuffA() {
    vec4 fragColor = vec4(0.0);

    vec4 userImage = vec4(0.);
	    vec4 clones = vec4(0.);
	    float colorKey = 0.;
	    float matte = 0.4-matte_tolerance*0.39;
	    vec3 matteColor = mix(color_key, texture(keyBuffer, vec2(0.0)).rgb, clamp(keying_style - 2., 0.0, 1.0));

		userImage = texture(BuffB, _uv);
		float mediaLuma = (dot(userImage.rgb, vec3(1.))/3.) + matte_tolerance;
        clones = userImage;
	    int numClones = int(iterations);
	    vec2 offset = vec2(0, 0.);
	    
	    float zoomMod = 1.5 - pow(clones_scale, 2.)*3. - 0.5*(sin(syn_BassTime*0.1) + cos(syn_HighTime*0.11 - 300. + TIME*0.05))*audio_reactive;
	    vec2 uv = _uv;
	    

        vec2 offsetMod = vec2(0.);
    	for(int i=0; i<numClones; i++) {
    	    uv = uv+(uv-0.5)*(zoomMod*i*.01);
    	    
    	    offsetMod = vec2(sin(syn_HighTime*0.05*PI*i), cos(syn_MidHighTime*0.03*PI*i))*0.07*audio_reactive;
    	    //xy
    	    offsetMod += sign(clones_offset)*vec2(pow(clones_offset, vec2(2.))*i*0.10);
    	    
    	    offset = vec2(uv) + offsetMod; // diaganol
    	    
    	   
    	    //kaleido
    	    float rep = 2.0*round(3.0 + syn_FadeInOut);
            float ss = (0.05*6.0 + dot(offsetMod*0.33*i, vec2(1.)))/rep;
            ss += NSIN(2.*i);
    	    vec2 kaleidoOffset = smoothKaleidoscope(vec2((_uvc)*sin(i*3.) - _uvc*sin(i*offsetMod*0.63)*0.3 + _uvc*cos(i*offsetMod*0.5) + (_uvc)*(zoomMod*i*.8)), ss, rep);
    	    offset = mix(offset, kaleidoOffset, clamp( cloner_style, 0.0, 1.0));
    	    
    	    //tunnel
    	    vec2 tunnelOffset = _toPolar(vec2(uv*2. - 1.)*(zoomMod*0.05*i) + offsetMod*i*0.03); // tunnel
    	    tunnelOffset = _rotate(tunnelOffset, dot(offsetMod, vec2(1.))*0.03*i);
    	    
    	    offset = mix(offset, tunnelOffset, clamp( cloner_style - 1.0, 0.0, 1.0));
    	    
    	    
    	    //turbulent noise
    	    vec2 turbulentNoiseOffset = vec2(uv - _fbm(sin(uv*3.*_uvc + syn_Time*0.1 + i) + offsetMod*0.3) + 0.5*uv*_noise(uv*3. + i*0.5 + 50.  + offsetMod) + offsetMod) +0.25;
    	    
    	    offset = mix(offset, turbulentNoiseOffset, clamp( cloner_style - 2.0, 0.0, 1.0));
    	    
    	    
       	    if(0.5 < edge_behavior && edge_behavior < 1.5) {
    	        offset.x = loop(offset.x);
    	        offset.y = loop(offset.y);
    	    }
    	    
    	    offset = mix(offset, clamp(offset, vec2(0.0), vec2(1.0)), clamp(edge_behavior - 1., 0.0, 1.0));
	    
    	   vec4 offsetClone = texture(BuffB, offset);
    	    clones = mix(offsetClone, clones, clones.a);
    	}
    	
    	
    //mix the clones with the base image
	fragColor = mix(clones, userImage, userImage.a);

    //show user media on click
    fragColor = mix(fragColor, _loadUserImage(), _click.x*clamp(keying_style - 2., 0.0, 1.0));


	return fragColor;
}

vec4 renderBuffB() {
    vec4 fragColor = vec4(0.);
    
        fragColor = texture(defaultImage, _uv);
        
        if(_exists(syn_UserImage)) {
            fragColor = _loadUserImage();
        }
        
        //vignette
        float vignette = getVignette();
        vec4 mask = getContentMask(fragColor);
        fragColor = mask*vignette;
        
        fragColor.a = mask.a;

    return fragColor;
}

vec4 renderKeyBuff() {
    vec4 fragColor = vec4(0.);
    
    vec4 selectedCol = texture(keyBuffer, vec2(0.));
    if(_click.x > 0.5) {
      selectedCol = _textureUserImage(_correctUserImageCoords(abs(_mouse.xy)/RENDERSIZE));
    }
    fragColor = selectedCol;
    return fragColor;
}

 vec4 renderOpticalFlowBuffer() {
     vec4 fragColor = vec4(0.);
     
     vec4 rawuserMedia = texture(BuffB, _uv);
     
     vec4 mediaBuildup = texture(opticalFlowBuffer, _uv);
     
     mediaBuildup +=  pow(rawuserMedia-0.1, vec4(2.))*0.1*matte_tolerance*smoothstep(0.04, 0.09, dot(rawuserMedia - 0.1, vec4(1.))/3.);
     mediaBuildup = mediaBuildup*0.993*smoothstep(0.04, 0.09, dot(rawuserMedia - 0.1, vec4(1.))/3.);
     
     fragColor = clamp(mediaBuildup, 0.0, 1.0);
     return fragColor;
 }



vec4 renderMain(void)
{
	
//media buffer
if (PASSINDEX == 0) {
 return renderBuffB();
}
//buffer to store user click color value
else if (PASSINDEX == 1) {
    return renderKeyBuff();
}
//feedback build up
else if (PASSINDEX == 2) {
    return renderOpticalFlowBuffer();
}
//clones buffer
else if (PASSINDEX == 3) {
    return renderBuffA();
}

}
