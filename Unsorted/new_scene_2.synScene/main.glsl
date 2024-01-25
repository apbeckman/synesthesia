vec4 renderMain(void)
{
	vec2 position = _uv;
	return vec4(position.x, slider1, position.y, 1.0);
}
