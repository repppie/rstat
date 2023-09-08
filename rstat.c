/* Heavily inspired by ministat */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <err.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

struct dataset {
	int n;
	char *name;
	double *vals;
	int vs;
	double mean;
	double ss;
};

int col = 1;
double conf = 0.95;
#define	MAX_DS 8
char sym[MAX_DS] = " x+*%#@O";
int nds;

double
integ(double lower, double upper, int m, double a,
    double (*f)(double x, double n))
{
	double h, s0, s1, s2;
	int i;

#if 1
	/* Simpson's rule */
	s1 = s2 = 0;
	s0 = (*f)(lower, a) + (*f)(upper, a);
	h = (upper - lower) / (2 * m);
	for (i = 1; i <= m * 2 - 1; i += 2)
		s1 += (*f)(lower + i * h, a);
	for (i = 2; i <= m * 2 - 2; i += 2)
		s2 += (*f)(lower + i * h, a);
	return (h / 3.0 * (s0 + 4.0 * s1 + 2.0 * s2));
#else
	h = (upper - lower) / m;
	s0 = 0;
	printf("h %lf\n", h);
#if 1
	/* Trapezoidal rule */
	for (i = 1; i < m; i++)
		s0 += (*f)(lower + i * h, a, b);
	return (h * (s0 + (*f)(lower, a, b) / 2.0 + (*f)(upper, a, b) / 2.0));
#else
	/* Midpoint rule */
	for (i = 0; i < m; i++)
		s0 += (*f)(lower + i * h + h / 2.0, a, b);
	return (h * s0);
#endif
#endif
}

/* t distribution pdf */
double
dt(double x, double n)
{
	return (tgamma((n+1.0)/2.0) / tgamma(n/2.0) *
	    pow(1.0+x*x/n, -(n+1.0)/2.0) / sqrt(M_PI * n));
}

/* t distribution cdf */
double
pt(double x, double dof)
{
	return (integ(-1000, x, 10000, dof, dt));
}

/* t distribution quantile */
double
qt(double alpha, double dof)
{
	double q;
	int i;

	q = 0;
	/* Newton's method solving pt(x, dof) = alpha */
	for (i = 0; i < 10; i++)
		q -= (pt(q, dof) - alpha) / dt(q, dof);
	return (q);
}

double
var(struct dataset *d)
{
	return (d->ss / (d->n - 1.0));
}

double
stddev(struct dataset *d)
{
	return (sqrt(var(d)));
}

double
max(double a, double b)
{
	return ((a > b) ? a : b);
}

double
min(double a, double b)
{
	return ((a < b) ? a : b);
}

double
quantile(struct dataset *d, double q)
{
	double frac, _i, m;
	int i, j;

	/* linear interpolation between points */
	m = (d->n - 1) * q;
	i = max(0, min(d->n - 1, floor(m)));
	j = max(0, min(d->n - 1, ceil(m)));
	frac = modf(m, &_i);

	return (d->vals[i] + (d->vals[j] - d->vals[i]) * frac);
}

int
cmp(const void *_a, const void *_b)
{
	double *a, *b;

	a = (void *)_a;
	b = (void *)_b;
	if (*a < *b)
		return (-1);
	else if (*a > *b)
		return (1);
	return (0);
}

struct dataset *
new_dataset(void)
{
	struct dataset *d;

	d = malloc(sizeof(struct dataset));
	memset(d, 0, sizeof(struct dataset));
	d->vs = 8;
	d->vals = malloc(d->vs * sizeof(double));
	return (d);
}

void
add_data(struct dataset *d, double v)
{
	double dt;

	if (d->n >= d->vs) {
		d->vs *= 2;
		d->vals = realloc(d->vals, d->vs * sizeof(double));
	}
	d->vals[d->n++] = v;
	/* Calculate mean and variance online */
	dt = v - d->mean;
	d->mean += dt / d->n;
	d->ss += dt * (v - d->mean);
}

struct dataset *
read_data(char *file)
{
	struct dataset *d;
	struct stat st;
	char *e, *end, *m, *s;
	int fd, i;

	if ((fd = open(file, O_RDONLY)) < 0)
		err(1, "open");
	if (fstat(fd, &st) != 0)
		err(1, "fstat");
	if ((m = mmap(NULL, st.st_size, PROT_READ | PROT_WRITE, MAP_PRIVATE, fd,
	    0)) == MAP_FAILED)
		err(1, "mmap");
	end = m + st.st_size;

	d = new_dataset();
	d->name = strdup(file);
	do {
		s = m;
		while (*m++ != '\n' && m < end);
		if (*s == '\n' || *s == '#')
			continue;
		*(m - 1) = '\0';
		i = 1;
		while ((e = strsep(&s, " \t")) != NULL) {
			if (*e == '\0')
				continue;
			if (i == col) {
				double v;
				v = strtod(e, NULL); /* XXX check error */
				add_data(d, v);
				break;
			}
			i++;
		}
	} while (m < end);

	qsort(d->vals, d->n, sizeof(double), cmp);
	printf("%c %s\n", sym[++nds], d->name);

	return (d);
}

/* Welch's t-test */
void
welch(struct dataset *d1, struct dataset *d2)
{
	double dof, p, q, se, se1, se2, t;

	se1 = sqrt(var(d1) / d1->n);
	se2 = sqrt(var(d2) / d2->n);
	se = sqrt(se1 * se1 + se2 * se2);
	t = d2->mean - d1->mean;
	t /= se;
	dof = pow(se, 4) / (pow(se1,4)/(d1->n-1) + pow(se2,4)/(d2->n-1));

	p = pt(-fabs(t), dof);
	/* Two tailed */
	if (2 * p > 1 - conf) {
		printf("No difference proven at 95%% confidence\n");
		exit(0);
	}

	q = qt(1 - (1 - conf)/2, dof);
	printf("diff: %g +/- %g [%g %g]\n", d2->mean - d1->mean, q * se,
	    (t - q) * se, (t + q) * se);
	/* XXX should be base be d1 or d2??? */
	printf("      %lf%% +/- %g%%\n", (d2->mean - d1->mean) * 100 / d1->mean,
	    q * se * 100 / d1->mean);
	printf("      (Welch's t %g p-val %g crit val %g se %g dof %g)\n", t,
	    2 * p, q, se, dof);
}

void
summary(struct dataset *d, char s)
{
	double iqr, q1, q3;
	int i, out, xout;

	q1 = quantile(d, 0.25);
	q3 = quantile(d, 0.75);
	iqr = q3 - q1;
	out = xout = 0;
	for (i = 0; i < d->n; i++) {
		if (d->vals[i] <= q1 - 3.0 * iqr || d->vals[i] >= q3 + 3.0 *
		    iqr)
			xout++;
		else if (d->vals[i] <= q1 - 1.5 * iqr || d->vals[i] >= q3 +
		    1.5 * iqr)
			out++;
	}
	printf("%c %3u %9.8g %9.8g %7.6g %7.6g %7.6g %7.6g %7.6g %4d"
	    " (%d mild)\n", s,
	    d->n, d->mean, stddev(d), quantile(d, 0), q1, quantile(d, 0.5), q3,
	    quantile(d, 1), xout, out);
}

int
main(int argc, char **argv)
{
	struct dataset *d1, *d2;

	if (argc < 3)
		errx(1, "Usage: %s <path> <path>\n", argv[0]);

	d1 = read_data(argv[1]);
	d2 = read_data(argv[2]);

	printf("    N      Mean    Stddev     Min     25p     50p     75p"
	    "     Max     Outliers\n");
	summary(d1, sym[1]);
	summary(d2, sym[2]);
	printf("\n");

	welch(d1, d2);

	return (0);
}
