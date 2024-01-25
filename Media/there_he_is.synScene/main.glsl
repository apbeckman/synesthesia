

			//******** BuffA Code Begins ********

const float spread = 20.0;
const float gravity = 0.15;
in vec3 textureDir; // direction vector representing a 3D texture coordinate
uniform samplerCube cubemap; // cubemap texture sampler
vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
    
    vec2 uv = fragCoord / RENDERSIZE.xy;
    vec2 pixel_uv = fragCoord;
    
    vec3 buffer = texture(BuffA, uv).rgb;
    
    vec3 off = vec3(1.0 / RENDERSIZE.x, 1.0 / RENDERSIZE.y, 0.0);
    
    vec3 normal = normalize(vec3(
        texture(BuffA, uv + off.xz).r - texture(BuffA, uv - off.xz).r,
        texture(BuffA, uv + off.zy).r - texture(BuffA, uv - off.zy).r,
        1.0));
    
    normal = normalize(normal + vec3(0.0, gravity, 0.0));
    
    float turbulence = pow(texture(image7, vec3(uv / 3.0 + TIME * 0.04, TIME * 0.05).xy).r, 1.25);
    float feedbackTurb = pow(texture(BuffA, vec3(uv / 3.0 + TIME * 0.04, TIME * 0.05).xy).r, 1.25);
    float turbMixed = mix(turbulence, feedbackTurb, 0.5);
    float val = texture(BuffA, uv + ((off.xy * normal.xy) * spread * turbulence)).r;
    
    val = val * 0.99;
    
    // input
    vec4 winput = texture(syn_UserImage, uv);
    float grey = (winput.r + winput.g + winput.b) / 3.0;
    grey = mix(grey, val, 0.125);
    winput.a *= step(0.5, distance(winput.rgb, vec3(0.0, 1.0, 0.0)));
    val = max(val, grey * winput.a);
    
    
    fragColor = vec4(vec3(clamp(val, 0.0, 1.0), grey, winput.a), 1.0);
    
	return fragColor; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord / RENDERSIZE.xy;
    vec3 c = texture(BuffA, uv).rgb;
    
    vec3 off = vec3(1.0 / RENDERSIZE.x, 1.0 / RENDERSIZE.y, 0.0);
    vec3 normal = normalize(vec3(
        texture(BuffA, uv + off.xz).r - texture(BuffA, uv - off.xz).r,
        texture(BuffA, uv + off.zy).r - texture(BuffA, uv - off.zy).r,
        c.r * c.r) + 0.005);
    
    
    vec3 col = texture(image7, uv + normal.xy).rgb * c.r;
    col += 0.65 * clamp(pow(dot(normal, vec3(0.0, -1.0, 0.75)), 10.0) * c.r, 0.0, 1.0);
    
    col = mix(col, vec3(c.g), c.b);
    
    
    fragColor = vec4(col, 1.0);
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