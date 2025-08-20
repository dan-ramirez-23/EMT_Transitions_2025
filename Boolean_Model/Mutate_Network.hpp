void Mutate_Network(RegNetwork *R)
{
	int edge1 = (int) ((R->numedges)*(distribution(generator)));
	int edge2 = (int) ((R->numedges)*(distribution(generator)));
	while(edge2 == edge1)
	{
		edge2 = (int) ((R->numedges)*(distribution(generator)));
	}
	int target1 = R->interactions[edge1][1];
	int target2 = R->interactions[edge2][1];

	R->interactions[edge1][1] = target2;
	R->interactions[edge2][1] = target1;

	return;
}
