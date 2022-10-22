mat2 rotate(float a) {
  return mat2(cos(a), -sin(a), sin(a), cos(a));
}

float smin( float a, float b, float k ) {
  float h = max( k-abs(a-b), 0.0 )/k;
  return min( a, b ) - h*h*h*k*(1.0/6.0);
}

float sphereSDF (vec3 p, float r) {
  return length(p) - r;
}

float displacement(vec3 p, float n) {
  return sin(p.x * n) * cos(p.y * n + n) - sin(p.z + p.y * n);
}

float map (vec3 p) {
  vec3 p1 = p;
  p1.xz *= rotate((script_time+terrain) * .75);
  p1.yz *= rotate((script_time+terrain) * .75);

  float s1 = sphereSDF(p1, .6);
  s1 += sin((p1.x + p1.y * p1.z) * (3.14 * 10.) - script_bass_time*0.245) * .015 - sin((p1.x - p1.y * p1.z) * (3.14 * 12.) - script_time*0.575) * .01;

  p.z -= (script_time*0.125)+rate;
  vec3 c = vec3(1.4, 0., 2.);
  p = mod(p, c) - .5 * c;

  float s2 = sphereSDF(p, .3);
  s2 += sin((p.x + p.y * p.z) * 20. + script_high_time*0.3) * .03 + cos(length(p.x - p.y * p.z) * 65. - script_time*0.46) * (.01+pow(syn_MidLevel*0.25, 2)*0.01);

  return smin(s1 * .5 + displacement(p1, 7. * (sin(p1.x * (p1.z + (smooth_hightime) * .125)) * .5 + .8)) * .1, s2 * .3, .3);
  //return smin(s1 * .5 + displacement(p1, 7. * (sin(p1.x * p1.z + iTime * .5) * .5 + .8)) * .1, s2 * .8, .3);
}

float trace (vec3 ro, vec3 rd) {
  float e = .001;
  float d = e * 2.;
  float t = 0.;
  for (int i = 0; i < 50; i++) {
    if (d < e || t > 60.) continue;
    d = map(ro + rd * t);
    t += d * (sin(TIME * .2) * .3 + .7);
    //t += d * (sin(iTime * .2) * .6 + .8);
  }
  return t;
}

vec3 getNormal (vec3 p) {
  float d = map(p);
  vec2 e = vec2(.01, 0.);

  return normalize(d - vec3(
    map(p - e.xyy),
    map(p - e.yxy),
    map(p - e.yyx)));
}

vec4 renderMainImage() 
{

	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

  vec2 uv = (fragCoord.xy / RENDERSIZE.xy) * 2.0 - 1.0;
  uv.x *= RENDERSIZE.x / RENDERSIZE.y;
  vec3 color;

  // camera
  vec3 ro = vec3(0., 1.+yPos, 3.-FOV);
  vec3 rd = normalize(vec3(uv.x, uv.y - 1., -3.));
  ro.xz *= rotate(TIME * .1);
  rd.xz *= rotate(TIME * .1);

  float t = trace(ro, rd);
  vec3 p = ro + rd * t;

  // texture
  color += sin(length(p + cos(atan(p.x, p.z * p.z) * 15.)) * 1.) * (t * .08);
  color *= cos(dot(p.xy, uv) * 5. - TIME) * .5;

  // lighting
  vec3 light = normalize(vec3(3., 2., 1.));
  light.xz *= rotate(smooth_hightime * .25)*0.3;
  vec3 nor = getNormal(p);
  float lOcclusion = .8;
  float ambient = clamp(.5 + .5 * nor.y, 0., 1.) * .3;
  float diffuse = clamp(dot(nor, light), 0., 1.) * .8;
  vec3 half_way = normalize(-rd + light);
  float specular = pow(clamp(dot(half_way, nor), 0.0, 1.0), 8.) * diffuse * 0.2;

  color += (ambient + diffuse + specular) * lOcclusion * vec3(.45, .6, .1);
  color += vec3(.13 + p.y * .5, .2 + p.x * .2, .8 + p.z * .2) * .6;

  float fog = 1. / (1. + t * t * .12);
  color *= fog * 2.;

  fragColor = vec4(color, 1.);
  return fragColor;
}

vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}