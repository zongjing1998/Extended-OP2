void saveMap(op_map m)
{
	char test[20];
	strcpy(test, "z_map_dat/");
	strcat(test, m->name);
	strcat(test, ".dat");
	FILE *fp;
	if ((fp = fopen(test, "w")) == NULL)
	{
		printf("can't open file %s\n", test);
		exit(2);
	}
	fprintf(fp, "name:%s\t from:%s\t from_size:%d\t to:%s\t to_size:%d\t dim:%d\n", m->name, m->from->name, m->from->size,
		m->to->name, m->to->size, m->dim);
	for (int i = 0; i < m->from->size; i++)
	{
		for (int j = 0; j < m->dim; j++)
		{
			fprintf(fp, "%d\t", m->map[i * m->dim + j]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}