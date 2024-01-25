
//Stack to apply all video aspect/metacontrol corrections
float smin(float a, float b, float k) {
    float h = max(k - abs(a-b), 0.) / k;
    return min(a, b) - h*h*h*k*1./6.;
}

vec4 noMediaPattern(vec2 uvIn) {

  vec2 p = 2.0*( _uv )-1.0;
  p.x *= RENDERSIZE.x/RENDERSIZE.y; 
  vec3 col = vec3(0);
  
  
  float c = sin(p.y*4.0+syn_CurvedTime+clamp(1.0/(30.0*abs(length(p.xy)-0.8)), 0.0, 1.0)*2.0);
  return vec4(c,(c>0.0)?c*c:0.0, 0.3, 1.0); 
}

/**
 * Detects edges using the Sobel equation
 * @name syn_pass_edgeDetectSobel
 * @param  {sampler2D} smp texture you wish to detect edges on
 * @returns {float} edges
 */

vec4 sobelIntensity(in vec4 color){
  return color;
  // return sqrt((color.x*color.x)+(color.y*color.y)+(color.z*color.z));
}
vec4 sobelHelper(float stepx, float stepy, vec2 center, sampler2D tex){
  // get samples around pixel
  vec4 tleft = sobelIntensity(texture(tex,clamp(center + vec2(-stepx,stepy), 0.0, 1.0)));
  vec4 left = sobelIntensity(texture(tex,clamp(center + vec2(-stepx,0), 0.0, 1.0)));
  vec4 bleft = sobelIntensity(texture(tex,clamp(center + vec2(-stepx,-stepy), 0.0, 1.0)));
  vec4 top = sobelIntensity(texture(tex,clamp(center + vec2(0,stepy), 0.0, 1.0)));
  vec4 bottom = sobelIntensity(texture(tex,clamp(center + vec2(0,-stepy), 0.0, 1.0)));
  vec4 tright = sobelIntensity(texture(tex,clamp(center + vec2(stepx,stepy), 0.0, 1.0)));
  vec4 right = sobelIntensity(texture(tex,clamp(center + vec2(stepx,0), 0.0, 1.0)));
  vec4 bright = sobelIntensity(texture(tex,clamp(center + vec2(stepx,-stepy), 0.0, 1.0)));

  vec4 x = tleft + 2.0*left + bleft - tright - 2.0*right - bright;
  vec4 y = -tleft - 2.0*top - tright + bleft + 2.0 * bottom + bright;
  vec4 color = sqrt((x*x) + (y*y));
  return color;
}
vec4 edgeDetectSobel(sampler2D tex){
	float stepSize = 1.0+expand*40.0;
  vec2 uv = _uv;
  if (uv.x < 0.0 || uv.y < 0.0 || uv.y > 1.0 || uv.x > 1.0) {
    return vec4(0.0);
  } 

  return sobelHelper(stepSize/RENDERSIZE.x, stepSize/RENDERSIZE.y, uv, tex);
}

/**
 * Shifts the RGB channels away from each other in the x and y axes.
 * @name syn_pass_rgbShift
 * @param  {sampler2D} smp texture you wish to affect
 * @param  {float} intensity how much you want to shift by.
 * @param  {float|vec2|vec3} speed speed of the shifted movement. 1 is a good default.
 * @returns {float} rgb shifted texture based on intensity
 */

vec4 rgbShift(sampler2D smp, float intensity, float speed){

  vec2 uv = _uv;

  float intensity_ = intensity*.1;
  float timeVariable = TIME*speed*.01;

  vec2 rPos = uv;
  vec2 gPos = uv;
  vec2 bPos = uv;

  if (woozy > 0.5) {
	  vec2 modifiedCenterR = vec2(
	    _statelessContinuousChaotic(timeVariable),
	    _statelessContinuousChaotic(timeVariable*1.3)
	    ) * intensity_;
	  vec2 modifiedCenterG = vec2(
	    _statelessContinuousChaotic(timeVariable*1.1),
	    _statelessContinuousChaotic(timeVariable*1.4)
	    ) * intensity_;
	  vec2 modifiedCenterB = vec2(
	    _statelessContinuousChaotic(timeVariable*1.2),
	    _statelessContinuousChaotic(timeVariable*1.5)
	    ) * intensity_;



	  rPos = (uv - modifiedCenterR*mix(1.0, syn_BassLevel*syn_BassLevel, auto_split));
	  gPos = (uv - modifiedCenterG*mix(1.0, syn_BassLevel*syn_BassLevel, auto_split));
	  bPos = (uv - modifiedCenterB*mix(1.0, syn_BassLevel*syn_BassLevel, auto_split));

	} else {
		rPos = (uv - _rotate(vec2(intensity_*(-1)*mix(1.0, syn_BassLevel*syn_BassLevel, auto_split), 0.0), syn_CurvedTime*0.1));
	  gPos = (uv - _rotate(vec2(intensity_*0*mix(1.0, syn_BassLevel*syn_BassLevel, auto_split), 0.0), syn_CurvedTime*0.1));
	  bPos = (uv - _rotate(vec2(intensity_*1*mix(1.0, syn_BassLevel*syn_BassLevel, auto_split), 0.0), syn_CurvedTime*0.1));
	}

  vec3 col;
  col.r = texture(smp, rPos).r;
  col.g = texture(smp, gPos).g;
  col.b = texture(smp, bPos).b;

  // vec2 minEdge = min(min(rPos, gPos), bPos);
  // vec2 maxEdge = max(max(rPos, gPos), bPos);

  // if (edges_fix > 0.5){
	 //  if (minEdge.x < 0.0 || minEdge.y < 0.0 || maxEdge.y > 1.0 || maxEdge.x > 1.0) {
	 //    col = vec3(0.0);
	 //  } 
  // }

 return vec4(col,1.0);
}


vec2 trianglePix(vec2 posIn, float sizeIn)
{
  vec2 tile_num = vec2(2.0,1.0)*sizeIn;
  vec2 uv = posIn;
  vec2 uv2 = floor(uv*tile_num)/tile_num;
  uv -= uv2;
  uv *= tile_num;
  vec2 posOut = uv2 + vec2(step(1.0-uv.y,uv.x)/(2.0*tile_num.x),                                          
                           step(uv.x,uv.y)/(2.0*tile_num.y));
  return posOut;
}

vec2 circlePix(vec2 posIn, float sizeIn)
{
        posIn *= RENDERSIZE;
        posIn.x *= 0.57735*2.0;
        posIn.y += mod(floor(posIn.x), 2.0)*0.5;
        posIn = abs((mod(posIn, 1.0) - 0.5));
        // posIn = vec2(posIn.x*1.5 + posIn.y, posIn.y*2.0);
        return posIn;
        // return abs(max(p.x*1.5 + p.y, p.y*2.0) - 1.0);
}
vec4 pixelate(sampler2D sampIn) 
{ 
  vec2 uv = _uv;

  float aspect = RENDERSIZE.y/RENDERSIZE.x;
  float sizeFixed = pow(pixel_size, 5.0);
  if (square_pix > 0.5){
    uv = vec2(_pixelate(uv.x, sizeFixed/aspect), _pixelate(uv.y, sizeFixed));
  } 
  if (tri_pix > 0.5){
    uv = trianglePix(_uv, sizeFixed);
  }
  vec4 col = vec4(0.0);
  if (circle_pix < 0.5){
    col = texture(sampIn, clamp(uv,0.0,1.0));
  } else {
    vec2 uv = _uv;
    vec2 imgSize = textureSize(sampIn, 0);
    float aspectImg = imgSize.x/imgSize.y;
    float sizeFixed = pow(pixel_size, 2.0)*0.8;
    vec2 smallGridSpace = mod(_uvc*4.5*pow(2, ceil(sizeFixed)), 2.0)-1.0;
    float smallXDiv = 4.0*pow(2, ceil(sizeFixed));
    float smallYDiv = 2.25*pow(2, ceil(sizeFixed));
    vec2 smallIndex = vec2(floor(uv.x*smallXDiv), floor((uv.y-0.5)*smallYDiv));
    // if (aspectImg < 16.0/9.0){
    //   smallXDiv *= aspect;
    // } else {
    //   smallYDiv *= aspect;
    // }
    // float smallXDivImg = aspectImg*9.0/4.0*pow(2, ceil(sizeFixed));
    // float smallYDivImg = 1/aspectImg*16.0/2.25*pow(2, ceil(sizeFixed));
    vec3 image = texture(sampIn, vec2(_pixelate(uv.x, smallXDiv), _pixelate(uv.y, smallYDiv))).rgb;
    float imgBri = dot(image, vec3(1.0))/3.0;

    float circSmall = 1.0-smoothstep((0.0+sqrt(imgBri*0.8*(0.5+syn_HighHits*syn_HighHits*0.7))), (0.05+sqrt(imgBri*0.8*(0.5+syn_HighHits*syn_HighHits*0.7))), 
    length(smallGridSpace));  
    col = vec4(vec3(circSmall*image), 1.0);
  }

  return col;
}
vec2 rotateCenter(vec2 uvIn, float amount){
  uvIn.y += (RENDERSIZE.x-RENDERSIZE.y)/RENDERSIZE.x;
  uvIn*=vec2(1.0, RENDERSIZE.y/RENDERSIZE.x);
  _uv2uvc(uvIn);
  uvIn = _rotate(uvIn, amount);
  _uvc2uv(uvIn);
  uvIn/=vec2(1.0, RENDERSIZE.y/RENDERSIZE.x);
  uvIn.y -= (RENDERSIZE.x-RENDERSIZE.y)/RENDERSIZE.x;
  return uvIn;
}

vec4 renderMain(void)
{
  if (PASSINDEX == 0){
		return _applyMediaContrast(_applyInvertMedia(noMediaPattern(_flipMediaCoords(_uv))));
  } else if (PASSINDEX == 1){
    vec4 img = vec4(0.0);
    if (!_isMediaActive()){
      img = texture(noMedia, _uv);
    } else {
      img = _loadMedia(); 
    }
    float mirror_x = smin(rd.x*-1, rd.x, 0.5);
    float mirror_y = smin(rd.y*-1, rd.y, 0.5);

    vec4 col = mix(img, img*extra_bright+texture(fxApplied, rotateCenter(_uv, fb_rotate*PI)-sign(fb_motion)*fb_motion*fb_motion*0.1+(_uvc)*fb_zoom*fb_zoom*sign(fb_zoom)*0.05), feedback_mix);
    return clamp_color > 0.5 ? clamp(col, 0.0, 1.0) : col;
  } else if (PASSINDEX == 2){
  	vec4 finalCol = vec4(0.0);

  	//Normal Image
  	vec4 normalCol = texture(userImgFB, _uv);
    vec4 edges = clamp(edgeDetectSobel(userImgFB), 0.0, 1.0);
    
    vec4 edgesCol = vec4(0.0);
    if (color_fix > 0.5){
      edgesCol = normalCol*edges*2.0;
    } else {
      edgesCol = edges;
    }

  	vec4 pixelateCol = pixelate(userImgFB);

  	finalCol += normalCol*(1.0-subtract*2.0)*normal_mix;

  	finalCol += edgesCol*edge_mix;

  	finalCol += pixelateCol*pixelate_mix;

    vec4 rgbCol = rgbShift(userImgFB, pow(split_amount,3.0)*0.5, woozy_speed);

    finalCol += rgbCol*rgb_mix;

  	return finalCol;
  } else if (PASSINDEX == 3.0){
    return texture(fxApplied, _uv);
  }
}
