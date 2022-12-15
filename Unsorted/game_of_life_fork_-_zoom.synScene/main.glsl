

			//******** Common Code Begins ********

float DigitBin(const in int x)
{
    return x==0?480599.0:x==1?139810.0:x==2?476951.0:x==3?476999.0:x==4?350020.0:x==5?464711.0:x==6?464727.0:x==7?476228.0:x==8?481111.0:x==9?481095.0:0.0;
}
float PrintValue(vec2 fragCoord, vec2 pixelCoord, vec2 fontSize, float value,
		float digits, float decimals) {
	vec2 charCoord = (fragCoord - pixelCoord) / fontSize;
	if(charCoord.y < 0.0 || charCoord.y >= 1.0) return 0.0;
	float bits = 0.0;
	float digitIndex1 = digits - floor(charCoord.x)+ 1.0;
	if(- digitIndex1 <= decimals) {
		float pow1 = pow(10.0, digitIndex1);
		float absValue = abs(value);
		float pivot = max(absValue, 1.5) * 10.0;
		if(pivot < pow1) {
			if(value < 0.0 && pivot >= pow1 * 0.1) bits = 1792.0;
		} else if(digitIndex1 == 0.0) {
			if(decimals > 0.0) bits = 2.0;
		} else {
			value = digitIndex1 < 0.0 ? fract(absValue) : absValue * 10.0;
			bits = DigitBin(int (mod(value / pow1, 10.0)));
		}
	}
	return floor(mod(bits / pow(2.0, floor(fract(charCoord.x) * 4.0) + floor(charCoord.y * 5.0) * 4.0), 2.0));
}
// Multiples of 4x5 work best
const vec2 fontSize = vec2(4,5) * vec2(5,5)*2.;
vec2 grid(int x, int y) {
    return fontSize.xx * vec2(1,ceil(fontSize.y/fontSize.x)) * vec2(x,y) + vec2(2);
}

#define mouseThreshold 5.

			//******** BuffA Code Begins ********

// load the texel color of a specific coordinate
vec3 loadValue(ivec2 p) {
    return texelFetch(BuffA, p, 0).rgb;
}

// count the alive neighbors surrounding a specific cell
int GetNeighbors(ivec2 p) {
    int count = 0;
    
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            if (i == 0 && j == 0) continue;
            if (loadValue(ivec2(p) + ivec2(i, j)).r > 0.5) {
                count ++;
            }
        }
    }
    
    return count;
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec3 color = vec3(0.);
    
    //bool isAlivePast = ( loadValue(ivec2(fragCoord)).r > 0.5 || 
        //ivec2(_mouse.xy) == ivec2(fragCoord) ) ? true : false;
    bool isAlivePast = ( loadValue(ivec2(fragCoord)).r > 0.5 )  ? true : false;
    int count = GetNeighbors(ivec2(fragCoord));
    
    // decide if the cell is alive or dead in this frame,
    // based on its status on the previous frame
    bool isAliveNow = false;
    if (isAlivePast && (count == 2 || count == 3)) {
        isAliveNow = true;
    } else if (!isAlivePast && count == 3) {
        isAliveNow = true;
    } else {
        isAliveNow = false;
    }
    
    if (isAliveNow) {
        color = vec3(1.);
    } else {
        color = vec3(0.);
    }
    
    if (_mouse.z > 0.5 && distance(_mouse.xy*0.5, fragCoord) < mouseThreshold) {
        color = vec3(2.);
    }
    
    fragColor = vec4(color, 1.);
	return fragColor; 
 } 


			//******** BuffB Code Begins ********

// load the texel color of a specific coordinate
vec3 loadChannel0(ivec2 p) {
    return texelFetch(BuffB, p, 0).rgb;
}

vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec3 color = loadChannel0(ivec2(fragCoord));
    if (_mouse.z > 0.5 && distance(_mouse.xy, fragCoord) < mouseThreshold) {
        color = vec3(1.);
    }
    
    fragColor = vec4(color, 1.0);
	return fragColor; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord/RENDERSIZE.xy*0.5;
    
    vec3 color = vec3(0.);
    //color = texture(BuffA, uv).xyz;
    color = texture(BuffA, uv).xyz;
    
    fragColor = vec4(color,1.);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderPassB();
	}
	if(PASSINDEX == 2){
		return renderMainImage();
	}
}