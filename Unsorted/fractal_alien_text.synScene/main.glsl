

			//******** BuffA Code Begins ********

/*

	draw letter shapes after subdividing uv space randomly

*/

#define PI 3.1415926535

float random2d(vec2 n) { 
    return fract(sin(dot(n, vec2(129.9898, 4.1414))) * 2398.5453);
}

vec2 getCellIJ(vec2 uv, float gridDims){
    return floor(uv * gridDims)/ gridDims;
}

vec2 rotate2D(vec2 position, float theta)
{
    mat2 m = mat2( cos(theta), -sin(theta), sin(theta), cos(theta) );
    return m * position;
}

//from https://github.com/keijiro/ShaderSketches/blob/master/Text.glsl
float letter(vec2 coord, float size)
{

    vec2 gp = floor(coord / size * 7.); // global
    
    vec2 rp = floor(fract(coord / size) * 7.); // repeated

    vec2 odd = fract(rp * 0.5) * 2.;
    float rnd = random2d(gp);
    float c = max(odd.x, odd.y) * step(0.5, rnd); // random lines
    c += min(odd.x, odd.y); // fill corner and center points
    c *= rp.x * (6. - rp.x); // cropping
    c *= rp.y * (6. - rp.y);
    return clamp(c, 0., 1.);
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;


    vec2 uv = fragCoord.xy / RENDERSIZE.xy;


    //correct aspect ratio
    uv.x *= RENDERSIZE.x/RENDERSIZE.y;


    float t = TIME;
    float scrollSpeed = 0.3;
    float dims = 2.0;
    int maxSubdivisions = 5;
    
    uv = rotate2D(uv,PI/12.0);
//    uv.y -= TIME * scrollSpeed;
    
    float cellRand;
    vec2 ij;
    
   	for(int i = 0; i <= maxSubdivisions; i++) { 
        ij = getCellIJ(uv, dims);
        cellRand = random2d(ij);
        dims *= 2.0;
        //decide whether to subdivide cells again
        float cellRand2 = random2d(ij + 454.4543);
        if (cellRand2 > 0.3){
        	break; 
        }
    }
   
    //draw letters    
    float b = letter(uv, 1.0 / (dims));
	
    //fade in
    float scrollPos = TIME*scrollSpeed + 0.5;
    float showPos = -ij.y + cellRand;
    float fade = smoothstep(showPos ,showPos + 0.05, scrollPos );
    b *= fade;
    
    //hide some
    //if (cellRand < 0.1) b = 0.0;
    
    fragColor = vec4(vec3(b), 1.0);
    
	return fragColor; 
 } 


/*

	Adding wobble and RGB shift PPO to Buffer A

*/

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord.xy / RENDERSIZE.xy;
   
    vec3 col = vec3(0.0);
    
    //wobble
    float strength = 0.008;
    float size = 5.0; 
    float speed = 0.3;
    vec2 p = -1.0 + 2.0 * uv;
    float time = TIME * speed;
    p.y -= time; //less wobble
    vec2 wobble = strength * vec2(cos(time+length(p*size)), sin(time+length(p*size)));
    uv += wobble;
    
    //offset rgb
    vec2 offset = vec2(2) / RENDERSIZE.xy;
    col.r = texture(BuffA,uv + offset).r;
    col.g = texture(BuffA,uv ).g;
    col.b = texture(BuffA,uv - offset).b;

    fragColor = vec4(col,1.0);
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