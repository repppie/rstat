/* Heavily inspired by ministat */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <err.h>
#include <time.h>
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
int boots, perms, seed;
double conf = 0.95;
#define	MAX_DS 8
char sym[MAX_DS] = " x+*%#@O";
int nds;

long
randint(int max)
{
	return ((max + 1) * (random() * 1.0 / RAND_MAX));
}

double
integ(double lower, double upper, int m, double a,
    double (*f)(double x, double n))
{
	double h, s0, s1, s2;
	int i;

	/* Simpson's rule */
	s1 = s2 = 0;
	s0 = (*f)(lower, a) + (*f)(upper, a);
	h = (upper - lower) / (2 * m);
	for (i = 1; i <= m * 2 - 1; i += 2)
		s1 += (*f)(lower + i * h, a);
	for (i = 2; i <= m * 2 - 2; i += 2)
		s2 += (*f)(lower + i * h, a);
	return (h / 3.0 * (s0 + 4.0 * s1 + 2.0 * s2));
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
free_dataset(struct dataset *d)
{
	/* XXX leaks name */
	free(d->vals);
	free(d);
}

struct dataset *
copy_dataset(struct dataset *old)
{
	struct dataset *d;

	d = malloc(sizeof(struct dataset));
	memcpy(d, old, sizeof(struct dataset));
	d->vals = malloc(d->vs * sizeof(double));
	memcpy(d->vals, old->vals, d->vs * sizeof(double));
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
	if (d->n < 3)
		errx(1, "%s needs at least 3 data points", d->name);

	return (d);
}

/* Destructively samples without replacement. */
struct dataset *
sample_no_repl(struct dataset *o, int n)
{
	struct dataset *d;
	double *v;
	int i, old_n, r;

	d = new_dataset();
	if (n >= o->n)
		errx(1, "sample_no_repl");
	for (i = 0; i < n; i++) {
		r = randint(o->n - 1);
		add_data(d, o->vals[r]);
		o->vals[r] = o->vals[--o->n];
	}
	v = malloc(o->vs * sizeof(double));
	memcpy(v, o->vals, o->n * sizeof(double));
	old_n = o->n;
	o->n = o->mean = o->ss = 0;
	for (i = 0; i < old_n; i++)
		add_data(o, v[i]);
	free(v);
	return (d);
}

/* Samples with replacement */
struct dataset *
sample_repl(struct dataset *o, int n)
{
	struct dataset *d;
	int i, r;

	d = new_dataset();
	for (i = 0; i < n; i++) {
		r = randint(o->n - 1);
		add_data(d, o->vals[r]);
	}
	return (d);
}

double
welch_tstat(struct dataset *d1, struct dataset *d2)
{
	double se, se1, se2, t;

	se1 = sqrt(var(d1) / d1->n);
	se2 = sqrt(var(d2) / d2->n);
	se = sqrt(se1 * se1 + se2 * se2);
	t = d2->mean - d1->mean;
	t /= se;
	return (t);
}

/*
 * Bootstrapped t-test without equality of variance.
 * See Algorithm 16.2 in Efron and Tibshirani (1993).
 */
void
bootstrap(struct dataset *d1, struct dataset *d2)
{
	struct dataset *b, *d, *dd1, *dd2, *n1, *n2;
	double diff, p, ql, qu, t;
	int i, g;

	t = welch_tstat(d1, d2);
	diff = d2->mean - d1->mean;
	d = copy_dataset(d1);
	for (i = 0; i < d2->n; i++)
		add_data(d, d2->vals[i]);
	n1 = new_dataset();
	for (i = 0; i < d1->n; i++)
		add_data(n1, d1->vals[i] - d1->mean + d->mean);
	n2 = new_dataset();
	for (i = 0; i < d2->n; i++)
		add_data(n2, d2->vals[i] - d2->mean + d->mean);
	b = new_dataset();

	g = 0;
	for (i = 0; i < boots; i++) {
		dd1 = sample_repl(n1, n1->n);
		dd2 = sample_repl(n2, n2->n);
		/* Two-tailed */
		if (fabs(welch_tstat(dd1, dd2)) >= fabs(t))
			g++;
		add_data(b, welch_tstat(dd1, dd2));
		free_dataset(dd1);
		free_dataset(dd2);
	}
	qsort(b->vals, b->n, sizeof(double), cmp);

	p = (g + 1.0) / (boots + 1.0);
	ql = quantile(b, (1 - conf) / 2);
	qu = quantile(b, 1 - (1 - conf) / 2);
	if (p > 1.0 - conf)
		printf("No difference proven at %.1f%% confidence\n", 100 *
		    conf);
	else {
		printf("Difference at %.1f%% confidence\n", 100 * conf);
		printf("      %g [%g %g]\n", diff, diff + ql * stddev(b),
		    diff + qu * stddev(b));
		printf("      %lf%% [%g%% %g%%]\n", diff * 100 / d1->mean,
		    (diff + ql * stddev(b)) * 100 / d1->mean,
		    (diff + qu * stddev(b)) * 100 / d1->mean);
	}
	printf("      (%d bootstrap samples, p-val %g t %g se %g seed %d)\n",
	    boots, p, t, stddev(b), seed);
}

/*
 * Permutation test
 * (In theory we don't gain anything by using the t-statistic instead of just
 * the difference in means. See problem 15.9 in Efron and Tibshirani (1993)).
 */
void
permute(struct dataset *d1, struct dataset *d2)
{
	struct dataset *d, *dd1, *dd2, *f;
	double diff, p, ql, qu, t;
	int i, g;

	diff = d2->mean - d1->mean;
	t = welch_tstat(d1, d2);
	d = copy_dataset(d1);
	for (i = 0; i < d2->n; i++)
		add_data(d, d2->vals[i]);
	f = new_dataset();

	g = 0;
	for (i = 0; i < perms; i++) {  
		dd2 = copy_dataset(d);
		dd1 = sample_no_repl(dd2, d1->n);
		/* Two-tailed */
		if (fabs(welch_tstat(dd1, dd2)) >= fabs(t))
			g++;
		add_data(f, welch_tstat(dd1, dd2));
		free_dataset(dd1);
		free_dataset(dd2);
	}
	qsort(f->vals, f->n, sizeof(double), cmp);

	/*
	 * XXX We can just look if p < alpha (not p * 2 < alpha), without
	 * looking if we're in the extreme quantiles.
	 */
	ql = quantile(f, (1 - conf) / 2);
	qu = quantile(f, 1 - (1 - conf) / 2);
	p = (g + 1.0) / (perms + 1.0);
	if (t > ql && t < qu)
		printf("No difference proven at %.1f%% confidence\n", 100 *
		    conf);
	else {
		printf("Difference at %.1f%% confidence\n", 100 * conf);
		printf("      %g [%g %g]\n", diff, diff + ql * stddev(f),
		    diff + qu * stddev(f));
		printf("      %lf%% [%g%% %g%%]\n", diff * 100 / d1->mean,
		    (diff + ql * stddev(f)) * 100 / d1->mean,
		    (diff + qu * stddev(f)) * 100 / d1->mean);
	}
	printf("      (%d permutations, p-val %g t %g se %g seed %d)\n", perms,
	    p, t, stddev(f), seed);
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
	/* Two-tailed */
	if (2 * p > 1 - conf) {
		printf("No difference proven at %.1f%% confidence\n", 100 *
		    conf);
		return;
	}

	q = qt(1 - (1 - conf)/2, dof);
	printf("Difference at %.1f%% confidence\n", 100 * conf);
	printf("      %g +/- %g [%g %g]\n", d2->mean - d1->mean, q * se,
	    (t - q) * se, (t + q) * se);
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

void
usage(void)
{
	fprintf(stderr, "Usage: rstat [-b samples] [-c confidence] [-C column]"
	    " [-p permutations] [-s seed] <file> [file ...]\n");
	fprintf(stderr, "\t-b : Bootstrapped t-test with n samples\n");
	fprintf(stderr, "\t-C : Column number to extract (starts at 1)\n");
	fprintf(stderr, "\t-c : Confidence in percent (defaults to 95)\n");
	fprintf(stderr, "\t-n : Print summary statistics only, no test\n");
	fprintf(stderr, "\t-p : Permutation test with n permutations\n");
	fprintf(stderr, "\t-s : Random seed\n");
	
	exit(1);
}

int
main(int argc, char **argv)
{
	struct dataset *d1, *d2;
	char c;
	int i, no_test;

	seed = time(NULL);

	no_test = 0;
	while ((c = getopt(argc, argv, "b:c:C:np:s:")) != -1) {
		switch (c) {
		case 'b':
			boots = atoi(optarg);
			break;
		case 'c':
			conf = strtod(optarg, NULL) / 100.0;
			if (conf <= 0 || conf >= 1)
				errx(1, "confidence needs to be in (0, 100)");
			break;
		case 'C':
			col = atoi(optarg);
			break;
		case 'n':
			no_test = 1;
			break;
		case 'p':
			perms = atoi(optarg);
			break;
		case 's':
			seed = atoi(optarg);
			break;
		default:
			usage();
		}
	}
	argc -= optind;
	argv += optind;

	if (argc < 1)
		usage();
	if (boots && perms)
		errx(1, "Only one of -b and -p can be set");
	srandom(seed);

	for (i = 0; i < argc; i++)
		printf("%c %s\n", sym[i+1], argv[i]);

	d1 = read_data(argv[0]);
	printf("    N      Mean    Stddev     Min     25p     50p     75p"
	    "     Max     Outliers\n");
	summary(d1, sym[++nds]);
	argv++;
	while (--argc) {
		d2 = read_data(argv[0]);
		summary(d2, sym[++nds]);
		if (!no_test) {
			if (boots)
				bootstrap(d1, d2);
			else if (perms)
				permute(d1, d2);
			else
				welch(d1, d2);
		}
		argv++;
	}

	return (0);
}
