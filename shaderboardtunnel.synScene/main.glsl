vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** Common Code Begins ********


#define TIME (TIME + 2.5)
    #define T TIME
    
    #define SPEED 2.

mat2 rot(float x) {
    return mat2(cos(x), sin(x), -sin(x), cos(x));
}

float half_tile(vec2 p, float k) {
    p.x -= 0.5;
    vec2 a = abs(p);
    float r1 = max(abs(p.x + 0.5), a.y) - 0.25;
    float r2 = length(vec2(p.x + 0.5, a.y)) - 0.25;
    float r = mix(r1, r2, k);
    vec2 s = (vec2(p.x, a.y) + vec2(0.25, 0.0)) * rot(3.141592 * 0.25);
    r = min(r, max(abs(s.x - 0.1) - 0.1, abs(s.y) - 0.01));
    return r;
}

float inner_tile(vec2 p) {
    p.x -= 0.5;
    return 0.23 - max(abs(p.x + 0.5), abs(p.y));
}

float tile(vec2 p) {
    p += vec2(0.5, -0.5);
    vec2 a = abs(p);
    float d = max(max(a.x, a.y) - 0.5, 0.48 - max(a.x, a.y));
    p.x += 0.5;
    p.y += 0.5;
	d = min(d, half_tile(p, 0.0));
    p.x -= 1.0;
    p.y -= 1.0;
    p.y = -p.y;
    p.x = -p.x;
    d = min(d, half_tile(p, 1.0));
    return d;
}

float rep_tile(vec2 p) {
	float d = 1000.0;
    d = min(d, tile(p));
    p.x = - p.x;
    p.y += 1.0;
    d = min(d, tile(p));
    return d;
}

float map(vec3 p) {
	vec3 a = abs(p);
    float d = -max(a.x - 2.0, a.y - 1.0);
    float az = (fract(p.z / 2.0) - 0.5) * 2.0;
    d = min(d, max(abs(p.x - 1.0), abs(az - 0.25)) - 0.25);
    d = min(d, max(abs(p.x + 1.0), abs(az + 0.25)) - 0.25);
    float bz = (fract((p.z + 1.0) / 2.0) - 0.5) * 2.0;
    d = min(d, max(abs(abs(p.y) - 1.0), abs(bz)) - 0.25);
    return d;
}

vec3 normal(vec3 p)
{
	vec3 o = vec3(0.0001, 0.0, 0.0);
    return normalize(vec3(map(p+o.xyy) - map(p-o.xyy),
                          map(p+o.yxy) - map(p-o.yxy),
                          map(p+o.yyx) - map(p-o.yyx)));
}

float trace(vec3 o, vec3 r) {
    float t = 0.0;
    for (int i = 0; i < 64; ++i) {
		t += map(o + r * t) * 0.7;
    }
    return t;
}

float thetex(vec2 p, float t) {
    p *= 4.0;
    vec2 pa = p - vec2(floor(mod(t, 2.0)), floor(mod(t, 3.0)));
    pa *= rot(3.141592 * 0.5 * floor(mod(t, 5.0)));
    vec2 s = mod(floor((pa + vec2(0.0, 3.0)) / 3.0), 2.0) * 2.0 - 1.0;
    pa = (fract(pa / 3.0) - 0.5) * 3.0;
    pa *= s.yx;
	float d = rep_tile(pa);
    d = min(d, abs(pa.x) - 0.02);
    d = max(d, 0.0);
    d = 1.0 / (1.0 + d * d * 10000.0);
    return d;
}

float vectex(vec3 p, float t) {
    float r = max(thetex(p.xy, t), thetex(p.xz, t + 1.0));
    return max(r, thetex(p.yz, t + 2.0));
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord / RENDERSIZE.xy;
    uv = uv * 2.0 - 1.0;
    uv.x *= RENDERSIZE.x / RENDERSIZE.y;
    
    vec3 r = normalize(vec3(uv, 1.0 - dot(uv, uv) * 0.2));
    vec3 o = vec3(0.0, -0.1, 1.0 + TIME + sin(3.41592 * TIME) * 0.12);

    float rt = clamp(fract(TIME), 0.0, 1.0);
    float rtt = clamp(fract(TIME), 0.0, 1.0);
    r.yz *= rot(sin(1.0 + TIME + sin(3.41592 * TIME) * 0.12) * 0.12);

    float t = trace(o, r);
    vec3 w = o + r * t;
    vec3 n = normal(w);
    float aoc = map(w + n * 0.3);

    float gt = 0.0;

    vec3 k = vec3(vectex(w, gt)) * vec3(1.0, 1.0, 0.5);
    
    vec3 col = vec3(0.0);
    col += k * 3.0 + (1.0 - k) * vec3(0.0, 0.25, 0.25);
    col += vec3(1.0, 0.0, 0.0) * abs(n.y);
    col *= 1.0 / (1.0 + pow(t * 0.4, 4.0));
    col *= mix(0.1, 1.0, aoc);
    
    fragColor = vec4(sqrt(col), 1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}