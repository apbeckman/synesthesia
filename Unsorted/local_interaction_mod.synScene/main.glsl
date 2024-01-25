

			//******** BuffA Code Begins ********
float dist(vec2 p0, vec2 pf) {
  return sqrt((pf.x - p0.x) * (pf.x - p0.x) + (pf.y - p0.y) * (pf.y - p0.y));
}
float d2 = dist(RENDERSIZE.xy * 0.5, _xy.xy) * 0.004;
float punch = punch_in - punch_out;

#define A(COORD) texture(BuffA,(COORD)/RENDERSIZE.xy)
#define C(COORD) texture(BuffC,(COORD)/RENDERSIZE.xy)
float mediaHits = mix(1.0, syn_MidHits, media_hits);
vec4 mediaEdges = texture(media_pass_fx, _uv);
vec4 media = texture(media_pass, _uv);

vec4 mediaEdge = vec4(mix(mediaEdges, media, 1.0 - mediaEdgeMix * 0.35));
float growthFactor = 0.25 + 2 * normalize(pow((syn_BassLevel * 0.5) + (syn_MidLevel * 0.25) + (syn_Intensity * 0.25), 2.0));
//vec4 image = texture(syn_UserImage,(_xy)/RENDERSIZE.xy);
vec4 image = vec4(0.);
vec3 logoCol = vec3(0.0);
float cosnoise = _scale(_fbm(cos(TIME * 0.5)), -1.0, 1.0);
float sinoise = _scale(_fbm(sin(TIME * 0.5)), -1.0, 1.0);
bool img = _exists(syn_UserImage);
vec4 mediaPass() {
  vec4 media = vec4(0.);
  vec2 U = _xy;
  media = _loadMedia();
  // media = mix(media, media*(low), media_hits);
  return media;
}
vec4 mediaPassFX() {
  vec4 media_fx = _edgeDetectSobel((media_pass));
  return media_fx;
}
vec4 feedback() {
    vec4 fb = _edgeDetectSobel(syn_FinalPass, _uv);
    vec2 U = _xy;
    return fb;
}
vec4 renderPassA() {
  vec4 Q = vec4(0.0);
  vec2 U = (_xy);
  U -= vec2(cosnoise, sinoise);
  mediaEdge *= (mediaEdge + 0.2);

  float scaleModder = Zoom * (.50 + 2.5 * growthFactor);
  U -= _uvc * scaleModder;

  if(_exists(syn_UserImage)) {
    logoCol = mediaEdge.rgb * media_impact;
  }
  float cosnoise = _scale(_fbm(cos(TIME)), -1.0, 1.0);
  float sinoise = _scale(_fbm(sin(TIME)), -1.0, 1.0);
  U += _uvc * Stretch + vec2(cosnoise, sinoise) * Stretch;
  U += vec2(_fbm(sin(_uvc.xy * PI))) * (Warp + 2 * vec2(_warp * sin(TIME * 0.2), _warp * cos(TIME * 0.2))) * (0.75 + 0.75 * growthFactor);
  U -= .5 * RENDERSIZE.xy;
  //  U *= .99975+(growthFactor*0.00625*Zoom);
  //  U *=(1.0+ 0.001*low*Zoom);
  U += (fisheye - punch * (1.0 + low)) * _uvc * d2 * (.50 + 2.5 * growthFactor);

    //U = _rotate(U, spin);
  float spin = (spin_left - spin_right) * PI;
  float a = .003 * sin(.125 * smoothTimeC - 2. * length(U - 0.5 * RENDERSIZE.xy) / RENDERSIZE.y) * (Spin + spin);
  U *= mat2(cos(a), -sin(a), sin(a), cos(a));
  U += .5 * RENDERSIZE.xy;
  U.xy += Dir_XY.xy * growthFactor * 0.5;
  image = vec4(logoCol, 0.);
    // image *= MediaImpact;
  image = mix(image, Q, 0.8 - syn_Level * 0.25);

  Q = A(U);

    // Neighborhood :
  vec4 pX = A(U + vec2(1, 0)) + (image * (0.25 + (0.25 * growthFactor))) * 0.375;
  vec4 pY = A(U + vec2(0, 1)) + (image * (0.25 + (0.25 * growthFactor))) * 0.375;
  vec4 nX = A(U - vec2(1, 0)) - (image * (0.25 + (0.25 * growthFactor))) * 0.375;
  vec4 nY = A(U - vec2(0, 1)) - (image * (0.25 + (0.25 * growthFactor))) * 0.375;
  vec4 m = 0.25 * (pX + nX + pY + nY);

  vec3 pulse = vec3(0.);
  vec2 circleposition = _muvc;
  float circle = 0.;
  vec2 cc = _uvc;
  float t = TIME;
  for(float i = 0.; i < 41; i++) {
    circle = length(cc);

    pulse += smoothstep(0.9, 0.4, circle) * vec3(_pulse(circle, fract(t + 3.33 * i), 0.1));
    vec2 circle_position = _uvc - circleposition;

    float center_dist = length(circle_position - i * 0.2);

    center_dist = smoothstep(0, 0.05, center_dist);

    pulse += vec3(center_dist);

  }

  float b = mix(1., abs(Q.z), .78);
  if(img) {
    Q = mix(image, Q, .50 - 0.25 * syn_Intensity * mediaHits);
  }
  Q.xyz += (1. - b + 0.0125 * basshits) * (0.25 * vec3(pX.z - nX.z, pY.z - nY.z, -pX.x + nX.x - pY.y + nY.y) - Q.xyz) * (Diffusion + growthFactor * 0.125);

  Q = mix(Q, m, b);
  if(length(Q.xy) > 0. || Reset > 0.)
    Q.xy = normalize(Q.xy);

  if(FRAMECOUNT <= 1 || Reset > 0.)
    Q = sin(.01 * length(U - (0.5 * RENDERSIZE.xy + PI * 0.25 * _uvc)) * vec4(1, 2, 3, 4));

  if(_mouse.z > 0. && length(U - _uvc - _mouse.xy) < (RENDERSIZE.y) * 0.125)
    Q *= 0.;
  //  Q.rgb += _fbm(Q.rgb)*0.001;

  return Q;
} 
/*bool Media = bool(media);
#ifdef Media
#define A1(COORD) texture(syn_UserImage,(COORD)/RENDERSIZE.xy)
#else*/
vec4 renderPassB() {
  // if(_exists(syn_UserImage)) {
  //   logoCol = mediaEdge.rgb * media_impact;
  // }

  vec4 Q = vec4(0.0);
  vec2 U = _xy;
  Q = texture(BuffA, _uv);
  // Q = mix(Q, Q + mediaEdge*media_impact, media_color_mix);
  Q = mix(Q, Q + media*media_impact, media_color_mix);

 // Q.rgb = mix(vec3(_fbm(Q.rgb)), Q.rgb, 1-0.01*syn_BassLevel);
  //Q = mix(image, Q, 0.9-syn_Hits*0.23);

//  Q += image*0.2;
  return Q;
}
vec4 renderPassC() {

  vec4 Q = vec4(0.0);
  vec2 U = _xy;
  Q = texture(BuffB, _uv);
  Q = mix((_edgeDetectSobel(BuffB, _uv)), Q, 1.0 - EdgeMix * 0.95);
  // Q = mix(Q, Q + mediaEdge, media_color_mix);
  Q = mix(Q, Q + media*media_impact, media_color_mix);

  return Q;
}

#define A1(COORD) texture(BuffB,(COORD)/RENDERSIZE.xy)
#define A2(COORD) texture(BuffB,(COORD)/RENDERSIZE.xy)
//  #endif
float ln(vec3 p, vec3 a, vec3 b) {
  return length(p - a - (b - a) * min(dot(p - a, b - a), 0.) / dot(b - a, b - a));
}
vec4 renderMainImage() {
  vec4 Q = vec4(0.0);
  vec2 U = _xy;
    //U-= _uvc*Zoom*(2*growthFactor);
  vec3 logoCol = vec3(0.0);
  if(_exists(syn_UserImage)) {
    logoCol = mediaEdge.rgb * (0.35 + 0.5 * syn_Intensity * mediaHits * media_impact);
  }
  Q = A1(U);

  // Q += _loadUserImageAsMask();
    //U+= _uvc*PI*0.1;
    // image *= MediaImpact;

    // vec4 pX  =  A1(U + vec2(1,0));
    // vec4 pY  =  A1(U + vec2(0,1));
    // vec4 nX  =  A1(U - vec2(1,0));
    // vec4 nY  =  A1(U - vec2(0,1));
  vec4 pX = A2(U + vec2(1, 0)) + (image * (0.25 + (0.25 * growthFactor))) * 0.375;
  vec4 pY = A2(U + vec2(0, 1)) + (image * (0.25 + (0.25 * growthFactor))) * 0.375;
  vec4 nX = A2(U - vec2(1, 0)) - (image * (0.25 + (0.25 * growthFactor))) * 0.375;
  vec4 nY = A2(U - vec2(0, 1)) - (image * (0.25 + (0.25 * growthFactor))) * 0.375;
  if(_mouse.z > 0. && length(U - _uvc - _mouse.xy) < (RENDERSIZE.y) * 0.125)
    Q *= 0.;

  vec3 n = normalize(vec3(pX.z - nX.z, pY.z - nY.z, 1));
   /// n += sin(_loadUserImage())*0.1;
  n += logoCol.rgb * syn_MidLevel * 0.85 * media_impact;
  vec3 r = reflect(n, vec3(0, 0, -2));
  r -= logoCol.rgb * syn_Presence * media_impact;
  if(img) {
    // Q = mix(image, Q, .50-syn_Level*0.25*syn_Presence*media_impact);
    Q += mix(image, vec4(0), .50 -  0.5 * syn_Presence * media_impact);
  }
  Q = (0.55 + 0.5 * sin(smoothTimeC * 0.125 + atan(Q.x, Q.y) * vec4(1, 2, 3, 4)));
  float d = ln(vec3(.4, .4, 12) * RENDERSIZE.xyy + _uvc.xyy, vec3(U, 0), vec3(U, 0) + r) / RENDERSIZE.y;
  d -= length(mediaEdge)*media_impact*0.5;
  Q = _brightness(_contrast(Q, 1.1), 0.995);
  Q *= exp(-d * d) * .5 + .5 * exp(-3. * d * d) + (syn_Intensity * pow(syn_HighLevel * 0.2 + syn_MidHighLevel * 0.2 + syn_Intensity * 0.2, 2.)) * Flash;
  Q = mix(Q, Q*0.65 +0.65* media*media_impact, media_color_mix*0.5);
   // Q =mix(Q, _edgeDetectSobel(Q));
  Q = _grayscale(Q, greyscale);
  return Q;
}

vec4 renderMain() {
  if(PASSINDEX == 0) {
    return renderPassA();
  }
  if(PASSINDEX == 1) {
    return renderPassB();
  }
  if(PASSINDEX == 2) {
    return renderPassC();
  }
  if(PASSINDEX == 3) {
    return mediaPass();
  }
  if(PASSINDEX == 4) {
    return mediaPassFX();
  }
  if(PASSINDEX == 5) {
    return feedback();
  }
  if(PASSINDEX == 6) {
    return renderMainImage();
  }
}