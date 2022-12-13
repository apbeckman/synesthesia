#define MAX_ITER 40

float termThres = 10.0;

vec3 filter_() {
  vec2 delta = 1. / RENDERSIZE;
  
  vec3 val = texture(backbuffer, _uv).xyz;


  vec3 l = texture(backbuffer, _uv + vec2(0., delta.y)).xyz;
  vec3 r = texture(backbuffer, _uv - vec2(0., delta.y)).xyz;
  vec3 u = texture(backbuffer, _uv + vec2(delta.x, 0.)).xyz;
  vec3 d = texture(backbuffer, _uv - vec2(delta.x, 0.)).xyz;

  vec3 n = vec3(_rand(_uvc+smoothTime)) - 0.5;
  
  vec3 bloom = max(val, max(max(l, r), max( u, d)));
  // bloom = bloom  + l + r + u + d;
  // bloom /= 5.; // orlando;
  return bloom + n/9.;

}

mat3 rotationMatrix(vec3 axis, float angle) {
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat3(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c);
}

vec3 opRep( vec3 p, vec3 c )
{
    vec3 q = mod(p, c) - 0.5 * c;
    
    return q;
}

float sdCappedCylinder( vec3 p, vec2 h )
{
  vec2 d = abs(vec2(length(p.xz),p.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float opS( float d1, float d2 )
{
    return max(-d1,d2);
}

float hex(vec3 p, vec2 h) {
  vec3 q = abs(p);
  return max(q.z-h.y, max((q.x*0.866025+q.y*0.5),q.y) - h.x);
}
float sdTriPrism( vec3 p, vec2 h )
{
    vec3 q = abs(p);
    return max(q.z-h.y,max(q.x*0.866025+p.y*0.5,-p.y)-h.x*0.5);
}
vec2 polar(vec2 c) {
  return vec2(atan(c.y, c.x), length(c));
}
vec2 cart(vec2 p) {
  return p.y * vec2(cos(p.x), sin(p.x));
}



vec2 DE(vec3 pos) {
  pos = pos * rotationMatrix(normalize(vec3(0., 0., 1.)), spin + TIME/10. + squiggle * sin(pos.z));
  
  // pos -= vec3(0., sin(TIME), cos(TIME));
  pos = opRep(pos, vec3(xy_mod.x+separation, xy_mod.y+separation,separation));
  

  pos.xy = polar(pos.xy);
  pos.x += 5.;
  pos.xy = cart(pos.xy);

  float a = hex(pos, vec2(size * separation,  depth/2.));
  float b = hex(pos, vec2(size_hole * size * separation, depth));

   // b = sdTriPrism(pos - vec3(-1., 0., 0.), vec2(0.6, 1.));
  // float c = sdTriPrism(pos - vec3(0., 0., 0.), vec2(0.6, 1.));

  float d =  max(a,-b);

  return vec2(pos.z, d);
}

vec3 gradient(vec3 p, float t) {
      vec2 e = vec2(0., t);

      return normalize( 
        vec3(
          DE(p+e.yxx).y - DE(p-e.yxx).y,
          DE(p+e.xyx).y - DE(p-e.xyx).y,
          DE(p+e.xxy).y - DE(p-e.xxy).y
        )
      );
    }

vec3 palette( in float t)
{
vec3 a = vec3(0.600, 0.500, 0.500);
vec3 b = vec3(1.000, 0.500, 0.500);
vec3 c = vec3(01.500, 0.690, 0.750);
vec3 d = vec3(0.500, 0.500, 0.500);
  
    return a + b*cos( 6.28318*(c*t+d) );
}

vec3 raycast() {
  
  vec3 camera = vec3( 0., 0., ((smoothTime)) * (1 + ( 0.5 + ((separation*1.1)-0.5) ) ) + fly_in_out);

  vec2 screenPos = -1.0 + 2.0 * _uv;
      screenPos  += (_uvc*PI*FOV*screenPos)*PI;

  screenPos.x *= RENDERSIZE.x / RENDERSIZE.y;

  vec2 n = vec2(_rand(_uvc+(smoothTimeC*0.4))) - 0.5;
  // screenPos += n/100.;
  vec3 ray = normalize(vec3( screenPos, 1.0));
  float thresh = exp(-termThres);
  // raycasting parameter 
  float t  = 0.;
  vec3 point;
  int iter = 0;
  bool hit = false;
  vec2 dist;
  // ray stepping 
  for(int i = 0; i < MAX_ITER+int(iters); i++) {
    point = camera + ray * t;
     dist = DE(point);

    thresh = exp(-termThres) * exp(t/4.);
  
    if (abs(dist.y) < thresh ) {
      hit = true;
      break;
    }
              
    t += dist.y * lighting;
    iter ++;

  }
          
  float shade = dot(gradient(point, 0.01 ), -ray);
  float ao = 1. -  float(iter) / float(MAX_ITER+iters);
  
  vec3 color = vec3(0.);
  
  if ( hit ){
    if (greyscale > 0.5){
      color = vec3(1.0);
    } else {
      color = _palette(_triWave(point.z/5.,3), 0.95*vec3(1.445, 1.107, 0.460), vec3(0.760, 1.000, 0.470), vec3(0.870, 0.341, 0.894), vec3(0.440, 0.404, 0.916));
    }
  }

  color *= sqrt(ao)*(0.5+shade*0.5);
  vec3 mediaCol = vec3(0.0);
  if (syn_MediaType > 0.5){
    if (dot(color, vec3(1.0))<0.1){
      color = _loadUserImage().rgb*(0.5+shade*1.7)*(length(vec2(0.5))-length(_uvc));
    }   
  }

  return color;


}

vec4 renderMain() {
  vec4 out_FragColor = vec4(0.0);

  if (PASSINDEX == 0) {
    out_FragColor = vec4(raycast(), 1.0);
  } else if (PASSINDEX == 1) {
    // out_FragColor = texture(backbuffer, _uv);
    out_FragColor = vec4(filter_(), 1.0);
  }

  return out_FragColor;
}
