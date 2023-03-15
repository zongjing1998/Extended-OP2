inline int type_error(const grid *a, const char *type)
{
	(void)a;
	return (strcmp(type, "grid") && strcmp(type, "grid:soa"));
}

inline int type_error(const elem *a, const char *type)
{
	(void)a;
	return (strcmp(type, "elem") && strcmp(type, "elem:soa"));
}

inline int type_error(const uvar *a, const char *type)
{
	(void)a;
	return (strcmp(type, "uvar") && strcmp(type, "uvar:soa"));
}

inline int type_error(const side *a, const char *type)
{
	(void)a;
	return (strcmp(type, "side") && strcmp(type, "side:soa"));
}

inline int type_error(const node *a, const char *type)
{
	(void)a;
	return (strcmp(type, "node") && strcmp(type, "node:soa"));
}

inline int type_error(const cvbd2d *a, const char *type)
{
	(void)a;
	return (strcmp(type, "cvbd2d") && strcmp(type, "cvbd2d:soa"));
}

inline int type_error(const para *a, const char *type)
{
	(void)a;
	return (strcmp(type, "para") && strcmp(type, "para:soa"));
}

template <class T>
inline int type_error(const T *a, const char *type)
{
	return type_error(*a, type);
}

