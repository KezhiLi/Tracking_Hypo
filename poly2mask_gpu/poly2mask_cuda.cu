__global__ void poly2mask_cuda(
	int   * mask,
	int     nMaskPoints,
	int     nPolygonEdges,
	float * xs,
	float * ys,
	int     height)
{
	int idx = blockDim.x * (gridDim.x * blockIdx.y + blockIdx.x) + threadIdx.x;
	if (idx >= nMaskPoints || nPolygonEdges < 3) // At least 3 polyon points.
	{
		return;
	}

	int x = idx / height;
	int y = idx % height;
	
	float x0, y0, x1, y1;
	int wn = 0;
	
    for (int i = 0; i < nPolygonEdges; i++)
    {
		x0 = xs[i];
		y0 = ys[i];
		
		x1 = xs[i+1];
		y1 = ys[i+1];
		
		if (y0 <= y && y < y1)
		{
			if (((x1 - x0) * (y - y0) - (x - x0) * (y1 - y0)) > 0)
			{
				++wn;
			}
		}
		else if (y1 <= y && y < y0)
		{
			if (((x1 - x0) * (y - y0) - (x - x0) * (y1 - y0)) < 0)
			{
				--wn;
			}
		}
		
    }
	
	if (wn != 0)
	{
		mask[idx] = 1;
	}
}
