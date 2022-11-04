

vec3 iResolution = vec3(RENDERSIZE, 1.);
float iTime = TIME;

#define PI 3.141592
#define orbs 30.




#define orbSize 6.46
#define sides 4.






vec4 orb(vec2 uv, float s, vec2 p, vec3 color, float c) {
  return pow(vec4(s / length(uv + p) * color, 1.), vec4(c));
}

mat2 rotate(float a) {
  return mat2(cos(a), -sin(a), sin(a), cos(a));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
  vec2 uv = (2. * fragCoord - iResolution.xy) / iResolution.y;
  uv *= zoom;
  uv /= dot(uv, uv);
  uv *= rotate(rotation * smoothTime / 40.);
  for (float i = 0.; i < orbs; i++) {
    uv.x += sinMul * sin(uv.y * yMul + smoothTimeC * xSpeed) + cos(uv.y / yDivide - (smoothTimeC/5.));
    uv.y += cosMul * cos(uv.x * xMul - smoothTime * ySpeed) - sin(uv.x / xDivide - (smoothTime/5.));
    float t = i * PI / orbs * 2.;
    float x = radius * tan(t);
    float y = radius * cos(t + (smoothTimeC) / 15.);
    vec2 position = vec2(x, y);
    vec3 color = cos(.02 * uv.x + .02 * uv.y * vec3(-2, 0, -1) * PI * 2. / 3. + PI * (float(i) / colorShift)) * 0.5 + 0.5;
    fragColor += .65 - orb(uv, orbSize, position, 1. - color, contrast);
  }
  fragColor.a = 1.0;
}

vec4 renderMain() { 
 	vec4 out_FragColor = vec4(0.0);

    mainImage(out_FragColor, _xy.xy);

return out_FragColor; 
 } 
