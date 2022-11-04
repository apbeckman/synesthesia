vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


#define i_time TIME
#define i_resolution RENDERSIZE

vec2 scene, eps = vec2(.001035, -.00035);
float td; // travel dist
float tt; // global time
float attr;

mat2 rot2(float r) {
  return mat2(cos(r), sin(r), -sin(r), cos(r));
}

float box(vec3 p, vec3 b) {
  vec3 q = abs(p) - b;
  return max(max(q.x, q.y), q.z);
}

float boxes(vec3 p) {
  vec3 q = p;

  q *= 0.6;
  q.xz *= rot2(sin(p.y*0.5)*0.4);

  q.x = mod(q.x + smooth_basstime*.65, (3.0-moreSquigglies));

  float d = box(abs(q) - vec3(1.5, 10.0, .8), vec3(.3, 10.0, .3));
  return d;
}

vec2 map(vec3 p) {

  // reverse sphere
  attr = max(0.0, (6.0 - (length(p) - 6.0)));

  vec2 res = vec2(boxes(p - vec3(.0, .0, 8.)), 5.0);

  // background sphere
  vec2 o = vec2(length(p - vec3(.0, .0, 20.))-10.0 + sin((p.x+smooth_basstime*(-1.6))*5.0)*.1, 15.0);

  res = (res.x < o.x) ? res : o;

  p += vec3(0, -.5, -2);
  p *= 0.9;

  // A
  vec3 q = p;

  float ttt = tt - 1.5;
  for (int i=0; i<2; i++) {
    q = abs(q) - vec3(0.3*attr, 0.2, 3);
    q.xy *= rot2(attr*(0.1+sin(ttt)*.2));
    q.xz *= rot2(attr*0.2+sin(ttt+1.0)*.2);
  }

  q.z = mod(q.z, 2.5);
  o = vec2(box(q - vec3(.0, 4, 0.), vec3(.5, .5, 10)), 10);

  res = (res.x < o.x) ? res : o;

  // B

  q = p;
  attr = max(0.0, (4.0 - (length(p) - 4.0)));

  ttt += .5;
  for (int i=0; i<2; i++) {
    q = abs(q) - vec3(0.3*attr, 0.2, 3);
    q.xy *= rot2(attr*(0.1+sin(ttt)*.3));
    q.xz *= rot2(attr*0.2+sin(ttt+1.0)*.3);
  }

  q.z = mod(q.z, 2.5);
  o = vec2(box(q - vec3(.0, 3, 0.), vec3(.5, .5, 10)), 10);

  res = (res.x < o.x) ? res : o;

  q *= 0.6;
  q.y -= 0.2;
  o = vec2(box(q - vec3(.1, 2.5, .0), vec3(.6, .6, 10)), 10);

  res = (res.x < o.x) ? res : o;

  q *= 0.6;
  q.y -= 0.2;
  o = vec2(box(q - vec3(.1, 2, .0), vec3(.7, .7, 10)), 10);

  res = (res.x < o.x) ? res : o;

  return res;
}

vec2 trace(vec3 ro, vec3 rd) {
  vec2 h, t = vec2(0.1);
  for (int i=0; i<128; i++) {
    h = map(ro+rd*t.x);
    if (h.x < .0001 || t.x > 100.0) break;
    t.x += h.x*0.6;
    t.y = h.y;
  }
  if (t.x > 100.0) t.x = 0.0;
  return t;
}

vec4 renderMainImage() {
	vec4 frag_color = vec4(0.0);
	vec2 frag_coord = _xy;

  vec2 uv = vec2(frag_coord.x / i_resolution.x, frag_coord.y / i_resolution.y);
  uv -= 0.5;
  uv /= vec2(i_resolution.y / i_resolution.x, 1);

  tt = mod(script_time*.5, 100.0);

  // ===============

  vec3 ro = vec3(0, 2, -10);
  vec3 rd = normalize(vec3(uv,0.5));

  vec3 ld = normalize(vec3(.3, .5, -.5));

  vec3 fog = vec3(0.25) * (1.0 - (length(uv)-0.2));
  vec3 col = fog;

  scene = trace(ro, rd);

  td = scene.x;
  if (td > 0.0) {
    vec3 pp = ro+rd*td;
    vec3 nor = normalize(eps.xyy*map(pp+eps.xyy).x +
                         eps.yyx*map(pp+eps.yyx).x +
                         eps.xyx*map(pp+eps.xyx).x +
                         eps.xxx*map(pp+eps.xxx).x);

    vec3 b0 = vec3(.1, .5, .8);
    vec3 b1 = vec3(0.7, .8, .8);
    vec3 base = mix(b0, b1, attr);

    float foggy;

    if (scene.y == 5.0) {
      foggy = -.00001*td*td*td;
    }

    if (scene.y == 10.0) {
      foggy = -.0000005*td*td*td;
    }

    if (scene.y == 15.0) {
      foggy = -.00005*td*td*td;
      base = vec3(0, 1, 0);
    }

    float diff = max(0.0, dot(nor, ld));

    float aor = td/30.;
    float ao = exp2(-2.*pow(max(0.,1.-map(pp+nor*aor).x/aor),2.));

    float fr = pow(1.+dot(nor,rd),2.);

    vec3 sss = vec3(0.5)*smoothstep(0.,1.,map(pp+ld*0.4).x/0.4);

    col = base*(0.9*ao+0.2)*(diff+sss);
    col = mix(col, b0, min(fr, 0.2));

    col = mix(col, fog, 1.0-exp(foggy));
  }


  col = pow(col, vec3(1.8)); // stronger colours
  col = vec3(1.0) - exp(-col * 3.0); // exposure hdr
  col = pow(col, vec3(.45)); // gamma

  frag_color = vec4(col, 1);
	return frag_color; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}