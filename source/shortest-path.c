/*
 * This file is part of an experimental software implementation for solving the
 * single-source shortest path problem in edge-weighted graphs. The algorithm
 * considered for the implementation is Dijkstra's algorithm for the shortest
 * path problem and it runs in linear time with respect to the size of the host
 * graph.
 *
 * The source code is configured for a gcc build for Intel microarchitectures.
 * Other builds are possible but it might require manual configuration of the
 * 'Makefile'.
 * 
 * The source code is subject to the following license.
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2018 Suhas Thejaswi
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#include<assert.h>
#include<sys/utsname.h>
#include<math.h>
#include<ctype.h>

/************************************************************* Configuration. */
#define BIN_HEAP
#define TRACK_OPTIMAL

#ifdef TRACK_RESOURCES
#include<omp.h>
#define TRACK_MEMORY
#define TRACK_BANDWIDTH
#endif

typedef int index_t; // default to 32-bit indexing

/********************************************************** Global constants. */

#define MAX_DISTANCE ((index_t)0x3FFFFFFF)
#define MATH_INF ((index_t)0x3FFFFFFF)
#define UNDEFINED -2

/************************************************************* Common macros. */

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

// Linked list navigation macros. 

#define pnlinknext(to,el) { (el)->next = (to)->next; (el)->prev = (to); (to)->next->prev = (el); (to)->next = (el); }
#define pnlinkprev(to,el) { (el)->prev = (to)->prev; (el)->next = (to); (to)->prev->next = (el); (to)->prev = (el); }
#define pnunlink(el) { (el)->next->prev = (el)->prev; (el)->prev->next = (el)->next; }
#define pnrelink(el) { (el)->next->prev = (el); (el)->prev->next = (el); }


/*********************************************************** Error reporting. */

#define ERROR(...) error(__FILE__,__LINE__,__func__,__VA_ARGS__);

static void error(const char *fn, int line, const char *func, 
                  const char *format, ...) 
{
    va_list args;
    va_start(args, format);
    fprintf(stderr, 
            "ERROR [file = %s, line = %d] "
            "%s: ",
            fn,
            (index_t)line,
            func);
    vfprintf(stderr, format, args);
    fprintf(stderr, "\n");
    va_end(args);
    abort();    
}

/********************************************************* Get the host name. */

#ifdef TRACK_RESOURCES
#define MAX_HOSTNAME 256

const char *sysdep_hostname(void)
{
    static char hn[MAX_HOSTNAME];

    struct utsname undata;
    uname(&undata);
    strcpy(hn, undata.nodename);
    return hn;
}
#endif

/********************************************** Memory allocation & tracking. */

#ifdef TRACK_MEMORY
#define MALLOC(x) malloc_wrapper(x)
#define CALLOC(x, y) calloc_wrapper((x), (y))
#define FREE(x) free_wrapper(x)

#else

#define MALLOC(x) malloc((x))
#define CALLOC(x, y) calloc((x), (y))
#define FREE(x) free((x))

#endif

#ifdef TRACK_MEMORY
index_t malloc_balance = 0;

struct malloc_track_struct
{
    void *p;
    size_t size;
    struct malloc_track_struct *prev;
    struct malloc_track_struct *next;
};

typedef struct malloc_track_struct malloc_track_t;

malloc_track_t malloc_track_root;
size_t malloc_total = 0;

#define MEMTRACK_STACK_CAPACITY 512
size_t memtrack_stack[MEMTRACK_STACK_CAPACITY];
index_t memtrack_stack_top = -1;

void *malloc_wrapper(size_t size)
{
    if(malloc_balance == 0) {
        malloc_track_root.prev = &malloc_track_root;
        malloc_track_root.next = &malloc_track_root;
    }
    void *p = malloc(size);
    if(p == NULL)
        ERROR("malloc fails");
    malloc_balance++;

    malloc_track_t *t = (malloc_track_t *) malloc(sizeof(malloc_track_t));
    t->p = p;
    t->size = size;
    pnlinkprev(&malloc_track_root, t);
    malloc_total += size;
    for(index_t i = 0; i <= memtrack_stack_top; i++)
        if(memtrack_stack[i] < malloc_total)
            memtrack_stack[i] = malloc_total;
    return p;
}

void *calloc_wrapper(size_t n, size_t size)
{
    if(malloc_balance == 0) {
        malloc_track_root.prev = &malloc_track_root;
        malloc_track_root.next = &malloc_track_root;
    }
    void *p = calloc(n, size);
    if(p == NULL)
        ERROR("malloc fails");
    malloc_balance++;

    malloc_track_t *t = (malloc_track_t *) malloc(sizeof(malloc_track_t));
    t->p = p;
    t->size = (n*size);
    pnlinkprev(&malloc_track_root, t);
    malloc_total += (n*size);
    for(index_t i = 0; i <= memtrack_stack_top; i++)
        if(memtrack_stack[i] < malloc_total)
            memtrack_stack[i] = malloc_total;
    return p;
}

void free_wrapper(void *p)
{
    malloc_track_t *t = malloc_track_root.next;
    for(;
        t != &malloc_track_root;
        t = t->next) {
        if(t->p == p)
            break;
    }
    if(t == &malloc_track_root)
        ERROR("FREE issued on a non-tracked pointer %p", p);
    malloc_total -= t->size;
    pnunlink(t);
    free(t);

    free(p);
    malloc_balance--;
}

index_t *alloc_idxtab(index_t n)
{
    index_t *t = (index_t *) MALLOC(sizeof(index_t)*n);
    return t;
}

void push_memtrack(void)
{
    assert(memtrack_stack_top + 1 < MEMTRACK_STACK_CAPACITY);
    memtrack_stack[++memtrack_stack_top] = malloc_total;
}

size_t pop_memtrack(void)
{
    assert(memtrack_stack_top >= 0);
    return memtrack_stack[memtrack_stack_top--];
}

size_t current_mem(void)
{
    return malloc_total;
}

double inGiB(size_t s)
{
    return (double) s / (1 << 30);
}

void print_current_mem(void)
{
    fprintf(stdout, "{curr: %.2lfGiB}", inGiB(current_mem()));
    fflush(stdout);
}

void print_pop_memtrack(void)
{
    fprintf(stdout, "{peak: %.2lfGiB}", inGiB(pop_memtrack()));
    fflush(stdout);
}

void inc_malloc_total(size_t size)
{
    malloc_total += size;
}

void dec_malloc_total(size_t size)
{
    malloc_total -= size;
}
#endif

/******************************************************** Timing subroutines. */

#ifdef TRACK_RESOURCES
#define TIME_STACK_CAPACITY 256
double start_stack[TIME_STACK_CAPACITY];
index_t start_stack_top = -1;

void push_time(void)
{
    assert(start_stack_top + 1 < TIME_STACK_CAPACITY);
    start_stack[++start_stack_top] = omp_get_wtime();
}

double pop_time(void)
{
    double wstop = omp_get_wtime();
    assert(start_stack_top >= 0);
    double wstart = start_stack[start_stack_top--];
    return (double) (1000.0*(wstop-wstart));
}
#endif

/********************************************************* Utility functions. */

/******************************************************* String manipulation. */
char* strlower( char *s)
{
   char* t;

   for(t = s; *s != '\0'; s++)
      *s = (char)tolower(*s);

   return(t);
}

char* itoa(index_t val)
{
    index_t base = 10;
	static char buf[32] = {0};
	index_t i = 30;
	for(; val && i ; --i, val /= base)
		buf[i] = "0123456789abcdef"[val % base];
	return &buf[i+1];
}

/**************************************************************  prefix sum. */
index_t prefixsum(index_t n, index_t *a, index_t k)
{
    index_t run = 0;
    for(index_t u = 0; u < n; u++) {
        index_t tv = a[u];
        a[u] = run;
        run += tv + k;
    }
    return run;
}


/*************************************************** Graph build subroutines. */

typedef struct graph
{
    index_t n;
    index_t m;
    index_t num_edges;
    index_t edge_capacity;
    index_t *edges;
} graph_t;

static index_t *enlarge(index_t m, index_t m_was, index_t *was)
{
    assert(m >= 0 && m_was >= 0);

    index_t *a = (index_t *) MALLOC(sizeof(index_t)*m);
    index_t i;
    if(was != (void *) 0) { 
        for(i = 0; i < m_was; i++) {
            a[i] = was[i];
        }
        FREE(was);
    }    
    return a;
}

graph_t *graph_alloc()
{
    graph_t *g = (graph_t *) MALLOC(sizeof(graph_t));
    g->n = 0; 
    g->m = 0; 
    g->num_edges      = 0;
    g->edge_capacity  = 100;
    g->edges          = enlarge(3*g->edge_capacity, 0, (void *) 0);
    
    return g;
}

void graph_free(graph_t *g)
{
    if(g->edges != NULL)
        FREE(g->edges);
    FREE(g);
}

void graph_add_edge(graph_t *g, index_t u, index_t v, index_t w)
{
    assert(u >= 0 && v >= 0 && u < g->n && v < g->n);

    if(g->num_edges == g->edge_capacity) {
        g->edges = enlarge(6*g->edge_capacity, 3*g->edge_capacity, g->edges);
        g->edge_capacity *= 2;
    }

    assert(g->num_edges < g->edge_capacity);

    index_t *e = g->edges + 3*g->num_edges;
    g->num_edges++;
    e[0] = u;
    e[1] = v;
    e[2] = w;
}

#define MAX_LINE_SIZE 1024

graph_t * graph_load(FILE *in)
{
#ifdef TRACK_RESOURCES
    push_time();
    push_memtrack();
#endif

    char buf[MAX_LINE_SIZE];
    char in_line[MAX_LINE_SIZE];
    index_t n = 0;
    index_t m = 0;
    index_t u, v, w;
    graph_t *g = graph_alloc();

    while(fgets(in_line, MAX_LINE_SIZE, in) != NULL) {
        char *line = strlower(in_line);
        int c = line[0];
        char *tok;
        strcpy(buf, line);
        tok = strtok(buf, " ");
        switch(c) {
        case 'e':
            if(!strcmp(tok, "edges")) {
                sscanf(line, "edges %d", &m);
                g->m = m;
                break;
            } 
            else if(!strcmp(tok, "e")) {
                sscanf(line, "e %d %d %d", &u, &v, &w);
                graph_add_edge(g, u-1, v-1, w);
                break;
            }
            break;
        case 'n':
            if(!strcmp(tok, "nodes")) {
                sscanf(line, "nodes %d", &n);
                g->n = n;
            }
            break;
        default:
            break;
        }
    }
    assert(g->n != 0);
    assert(g->m == g->num_edges && g->m != 0);

#ifdef TRACK_RESOURCES
    double time = pop_time();
    fprintf(stdout, "input: n = %d, m = %d [%.2lfms] ",g->n, g->m, time);
    print_pop_memtrack();
    fprintf(stdout, " ");
    print_current_mem();
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    return g;
}

/******************************************************** Root query builder. */

typedef struct steinerq
{
    index_t  n;
    index_t  m;
    index_t  src;
    index_t  dst;
    index_t *pos;
    index_t *adj;
}dijkstra_t;

dijkstra_t *root_build(graph_t *g, index_t src, index_t dst)
{
#ifdef TRACK_RESOURCES
    fprintf(stdout, "root build: ");
    push_time();
    push_memtrack();
#endif
    index_t n    = g->n;
    index_t m    = g->m;
    index_t *pos = (index_t *) MALLOC((n+1)*sizeof(index_t));
    index_t *adj = (index_t *) MALLOC(((n+1)+(4*m)+(2*n))*sizeof(index_t));

    assert(src >= 0 && src < n);
    assert(dst >= -1 && dst < n);

    dijkstra_t *root = (dijkstra_t *) MALLOC(sizeof(dijkstra_t));
    root->n   = n;
    root->m   = m;
    root->src = src;
    root->dst = dst;
    root->pos = pos;
    root->adj = adj;

#ifdef TRACK_RESOURCES
    push_time();
#endif
    for(index_t u = 0; u < n; u++)
        pos[u] = 0;
#ifdef TRACK_RESOURCES
    double time = pop_time();
    fprintf(stdout, "[zero: %.2lfms] ", time);
    fflush(stdout);
    push_time();
#endif

    index_t *e = g->edges;
    for(index_t j = 0; j < 3*m; j+=3) {
        pos[e[j]]+=2;
        pos[e[j+1]]+=2;
    }

    pos[n] = (2*n);
    index_t run = prefixsum(n+1, pos, 1);
    assert(run == ((n+1)+(4*m)+(2*n)));

#ifdef TRACK_RESOURCES
    time = pop_time();
    fprintf(stdout, "[pos: %.2lfms] ", time);
    fflush(stdout);
    push_time();
#endif

    for(index_t u = 0; u < n; u++)
        adj[pos[u]] = 0;

    for(index_t j = 0; j < 3*m; j+=3) {
        index_t u    = e[j+0];
        index_t v    = e[j+1];
        index_t w    = e[j+2];
        index_t pu   = pos[u];
        index_t pv   = pos[v];
        index_t i_pu = pu + 1 + (2*adj[pu]);
        index_t i_pv = pv + 1 + (2*adj[pv]);

        adj[i_pv]     = u;
        adj[i_pv + 1] = w;
        adj[pv]++;
        adj[i_pu]     = v;
        adj[i_pu + 1] = w;
        adj[pu]++; 
    }

    index_t u = n;
    adj[pos[u]] = 0;
    index_t pu = pos[u];
    for(index_t v = 0; v < n; v++) {
        index_t i_pu  = pu + 1 + (2*adj[pu]);
        adj[i_pu]     = v;
        adj[i_pu + 1] = MATH_INF;
        adj[pu]++;
    }

#ifdef TRACK_RESOURCES
    time = pop_time();
    fprintf(stdout, "[adj: %.2lfms] ", time);
    time = pop_time();
    fprintf(stdout, "done. [%.2lfms] ", time);
    print_pop_memtrack();
    fprintf(stdout, " ");
    print_current_mem();
    fprintf(stdout, "\n");
    fflush(stdout);
#endif

    return root;
}

void dijkstra_free(dijkstra_t *root)
{
    if(root->pos != NULL)
        FREE(root->pos);
    if(root->adj != NULL)
        FREE(root->adj);
    FREE(root);
}

/********************************************************* Debug routines. */

#ifdef DEBUG
void print_adj(index_t u, index_t *pos, index_t *adj)
{
    index_t p = pos[u];
    index_t nu = adj[p];
    index_t *adj_u = adj + p + 1;
    fprintf(stdout, "adjacency list (%d) : ", u);
    for(index_t i = 0; i < nu; i++)
        fprintf(stdout, " %d %d|", adj_u[2*i]+1, adj_u[2*i+1]);
    fprintf(stdout, "\n");
}

void print_dist(index_t n, index_t *d)
{
    fprintf(stdout, "Shortest distance: \n");
    for(index_t u = 0; u < n; u++)
        fprintf(stdout, "%d: %d\n", u+1, d[u]);
    fflush(stdout);
}


void print_dist_matrix(index_t *d_N, index_t n)
{
    fprintf(stdout, "distance matrix: \n");
    for(index_t u = 0; u < n; u++) {
        for(index_t v = 0; v < n; v ++) {
            fprintf(stdout, " %d", d_N[u*n+v]);
        }
        fprintf(stdout, "\n");
    }
    fflush(stdout);
}

void print_graph(graph_t *g)
{
    fprintf(stdout, "graph_t: \n");
    fprintf(stdout, "n: %d\n", g->n);
    fprintf(stdout, "m: %d\n", g->m);
    fprintf(stdout, "num edges: %d\n", g->num_edges);
    fprintf(stdout, "edge capacity: %d\n", g->edge_capacity);
    fprintf(stdout, "edges: \n");
    for(index_t i = 0; i < g->num_edges; i++) {    
        index_t *e = g->edges + (3*i);
        index_t u = e[0];
        index_t v = e[1];
        index_t w = e[2];
        fprintf(stdout, "E %d %d %d\n", u+1, v+1, w);
    }
    fflush(stdout);
}

void print_dijkstra_t(dijkstra_t *root)
{
    fprintf(stdout, "steinerq_t: \n");
    fprintf(stdout, "n: %d\n", root->n);
    fprintf(stdout, "m: %d\n", root->m);
    index_t *pos = root->pos;
    index_t *adj = root->adj;
    fprintf(stdout, "pos:");
    for(index_t i = 0; i < root->n; i++)
        fprintf(stdout, " %d", pos[i]);
    fprintf(stdout, "\nadj:\n");
    index_t n = root->n + 1;
    for(index_t u = 0; u < n; u++) {
        index_t pu = pos[u];
        index_t adj_u = adj[pu];
        fprintf(stdout, "node: %d edges: %d|", u+1, adj_u);
        for(index_t i = 0; i < adj_u; i++) {
            fprintf(stdout, " %d %d|", 
                            adj[pu + 1 + (2*i)]+1, 
                            adj[pu + 1 + (2*i+1)]);
        }
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
    fflush(stdout);
}
#endif 

/******************************************************* Heap implementaions. */

/************************************************ Binary heap implementation. */
#ifdef BIN_HEAP
typedef struct bheap_item
{
    index_t item;
    index_t key; 
} bheap_item_t;

typedef struct bheap
{
    index_t max_n;
    index_t n;       // size of binary heap
    bheap_item_t *a; // stores (distance, vertex) pairs of the binary heap
    index_t *p;      // stores the positions of vertices in the binary heap
#ifdef TRACK_BANDWIDTH
    index_t key_comps;
    index_t mem;
#endif
} bheap_t;

bheap_t * bh_alloc(index_t n)
{
    bheap_t *h = (bheap_t *) malloc(sizeof(bheap_t));
    h->max_n = n;
    h->n = 0; 
    h->a = (bheap_item_t *) malloc((n+1)*sizeof(bheap_item_t));
    h->p = (index_t *) malloc(n*sizeof(index_t));
#ifdef TRACK_BANDWIDTH
    h->key_comps  = 0; 
    h->mem = 0;
#endif
#ifdef TRACK_MEMORY
    inc_malloc_total(sizeof(bheap_t) + 
                     ((n+1)*sizeof(bheap_item_t)) +
                     (n*sizeof(index_t)));
#endif
    return h;
}

void bh_free(bheap_t *h)
{
#ifdef TRACK_MEMORY
    index_t n = h->max_n;
    dec_malloc_total(sizeof(bheap_t) + 
                     ((n+1)*sizeof(bheap_item_t)) +
                     (n*sizeof(index_t)));
#endif
    free(h->a);
    free(h->p);
    free(h);
}

/************************************************** Binary heap operations. */

static void bh_siftup(bheap_t *h, index_t p, index_t q)
{
    index_t j = p;
    index_t k = 2 * p;
    bheap_item_t y = h->a[p];
#ifdef TRACK_BANDWIDTH
    index_t mem = 0;
    index_t key_comps = 0;
#endif

    while(k <= q) {
        bheap_item_t z = h->a[k];
        if(k < q) {
#ifdef TRACK_BANDWIDTH
            mem++;
            key_comps++;
#endif
            if(z.key > h->a[k + 1].key) z = h->a[++k];
        }

#ifdef TRACK_BANDWIDTH
        mem += 2;
        key_comps++;
#endif
        if(y.key <= z.key) break;
        h->a[j] = z;
        h->p[z.item] = j;
        j = k;
        k = 2 * j;
    }

    h->a[j] = y;
    h->p[y.item] = j;

#ifdef TRACK_BANDWIDTH
    h->mem += mem;
    h->key_comps += key_comps;
#endif
}

bheap_item_t bh_min(bheap_t *h)
{
    return (bheap_item_t) h->a[1];
}

static void bh_insert(bheap_t *h, index_t item, index_t key)
{
    index_t i = ++(h->n);
#ifdef TRACK_BANDWIDTH
    index_t mem = 0;
    index_t key_comps = 0;
#endif

    while(i >= 2) {
        index_t j = i / 2;
        bheap_item_t y = h->a[j];

#ifdef TRACK_BANDWIDTH
        mem ++;
        key_comps++;
#endif
        if(key >= y.key) break;

        h->a[i] = y;
        h->p[y.item] = i;
        i = j;
    }

    h->a[i].item = item;
    h->a[i].key = key;
    h->p[item] = i;
#ifdef TRACK_BANDWIDTH
    h->mem += mem;
    h->key_comps += key_comps;
#endif
}

static void bh_delete(bheap_t *h, index_t item)
{
    index_t n = --(h->n);
    index_t p = h->p[item];
#ifdef TRACK_BANDWIDTH
    index_t mem = 0;
    index_t key_comps = 0;
#endif

    if(p <= n) {
#ifdef TRACK_BANDWIDTH
        key_comps++;
        mem += 2;
#endif
        if(h->a[p].key <= h->a[n + 1].key) {
            h->a[p] = h->a[n + 1];
            h->p[h->a[p].item] = p;
            bh_siftup(h, p, n);
        }
        else {
            h->n = p - 1;
            bh_insert(h, h->a[n + 1].item, h->a[n+1].key);
            h->n = n;
        }
    }
#ifdef TRACK_BANDWIDTH
    h->mem += mem;
    h->key_comps += key_comps;
#endif
}

static void bh_decrease_key(bheap_t *h, index_t item, index_t new_key)
{
#ifdef TRACK_BANDWIDTH
    index_t mem = 1;
    index_t key_comps = 0;
#endif

    index_t i = h->p[item];
    while(i >= 2) {
        index_t j = i / 2;
        bheap_item_t y = h->a[j];

#ifdef TRACK_BANDWIDTH
        mem ++;
        key_comps++;
#endif
        if(new_key >= y.key) break;

        h->a[i] = y;
        h->p[y.item] = i;
        i = j;
    }

    h->a[i].item = item;
    h->a[i].key = new_key;
    h->p[item] = i;
#ifdef TRACK_BANDWIDTH
    h->mem += mem;
    h->key_comps += key_comps;
#endif
}

static index_t bh_delete_min(bheap_t * h)
{    
    bheap_item_t min = (bheap_item_t) h->a[1];
    index_t u = min.item;
    bh_delete((bheap_t *)h, u);
    return u;
}
#endif


/************************************************** Heap wrapper functions. */

#ifdef BIN_HEAP
// allocation
#define heap_alloc(n) bh_alloc((n))
#define heap_free(h) bh_free((bheap_t *)(h));
// heap operations
#define heap_insert(h, v, k) bh_insert((h), (v), (k))
#define heap_delete_min(h) bh_delete_min((h));
#define heap_decrease_key(h, v, k) bh_decrease_key((h), (v), (k));
// fetch structure elements
#define heap_n(h) ((bheap_t *)h)->n;
#define heap_key_comps(h) ((bheap_t *)h)->key_comps;
#define heap_mem(h) (h)->mem;
// heap nodes
#define heap_node_t bheap_item_t
#define heap_t bheap_t
#endif

/************************************************************ traceback path. */

void tracepath(index_t n, index_t s, index_t cost, 
               index_t v, index_t *p)
{
#ifdef TRACK_RESOURCES
    fprintf(stdout, "\ntracepath: ");
#endif
    fprintf(stdout, "[source: %d] [destination: %d] [cost: %s]\n", 
                     s+1, v+1, cost==MAX_DISTANCE?"INFINITY":itoa(cost));
    if(s == v) {
        fprintf(stdout, "\n\n");
        return;
    }

    if(cost == MAX_DISTANCE) {
        fprintf(stdout, "No path from %d to %d\n\n", s+1, v+1);
        return;
    }

    index_t u = p[v];
    while(u != s) {
        if(u == UNDEFINED)
            return;
        fprintf(stdout, "E %d %d\n", v+1, u+1);
        v = u;
        u = p[v];
    }
    fprintf(stdout, "E %d %d\n\n", v+1, u+1);
}

/*************************************************** Dijkstra shortest path. */

void dijkstra(index_t n,
              index_t m, 
              index_t *pos, 
              index_t *adj, 
              index_t s, 
              index_t *d,
              index_t *visit,
              index_t *p
#ifdef TRACK_BANDWIDTH
              ,index_t *heap_ops
#endif
             )
{
#ifdef TRACK_RESOURCES
    fprintf(stdout, "dijkstra: ");
    push_time();
    push_time();
#endif

    heap_t *h = heap_alloc(n);
    //initialise
    for(index_t v = 0; v < n; v++) {
        d[v]     = MAX_DISTANCE; //mem: n
        visit[v] = 0;            //mem: n
        p[v]     = UNDEFINED;    //mem: n
    }
    d[s] = 0;
    p[s] = UNDEFINED;

    for(index_t v = 0; v < n; v++)
        heap_insert(h, v, d[v]);



#ifdef TRACK_RESOURCES
    fprintf(stdout, "[init: %.2lfms]", pop_time());
    push_time();
#endif
    //visit and label
    while(h->n > 0) {
        index_t u = heap_delete_min(h); 
        visit[u]  = 1;

        index_t pos_u  = pos[u];
        index_t *adj_u = adj + pos_u;
        index_t n_u    = adj_u[0];
        for(index_t i = 1; i <= 2*n_u; i += 2) {
            index_t v   = adj_u[i];
            index_t d_v = d[u] + adj_u[i+1];
            if(!visit[v] && d[v] > d_v) {
                d[v] = d_v;
                p[v] = u;
                heap_decrease_key(h, v, d_v);
            }
        }
        //mem: 2n+6m
    }
#ifdef TRACK_RESOURCES
    fprintf(stdout, " [visit: %.2lfms] [%.2lfms]\n", pop_time(), pop_time());
    fflush(stdout);
#endif


#ifdef TRACK_BANDWIDTH
    *heap_ops = heap_mem(h);
#endif
    //free heap
    heap_free(h);
}


index_t dijkstra_query(dijkstra_t *root)
{
#ifdef TRACK_RESOURCES
    push_memtrack();
    push_time();
#endif

    index_t n        = root->n;
    index_t m        = root->m;
    index_t src      = root->src;
    index_t dst      = root->dst;
    index_t *d       = (index_t *) MALLOC(n*sizeof(index_t));
    index_t *visit   = (index_t *) MALLOC(n*sizeof(index_t));
    index_t *p       = (index_t *) MALLOC(n*sizeof(index_t));
#ifdef TRACK_BANDWIDTH
    index_t heap_ops = 0;
#endif

    //execute query
#ifdef TRACK_RESOURCES
    push_time();
#endif
    dijkstra(n, m, root->pos, root->adj, src, d, visit
            ,p
#ifdef TRACK_BANDWIDTH
            ,&heap_ops
#endif
            );

#ifdef TRACK_RESOURCES
    double dijkstra_time = pop_time();
    push_time();
#endif

    //trace path
    index_t min_cost = 0;
    if(dst != -1) {
        min_cost = d[dst];
        tracepath(n, src, min_cost, dst, p);
    } else {
        for(index_t v = 0; v < n; v++) {
            if(v != src) {
                min_cost = d[v];
                tracepath(n, src, min_cost, v, p);
            }
        }
    }

#ifdef TRACK_RESOURCES
    double trace_time = pop_time();
    index_t mem_graph = 5*n + 6*m;
    index_t trans_bytes = (mem_graph * sizeof(index_t) + 
                           heap_ops * sizeof(heap_node_t));
    double trans_rate  = trans_bytes / (dijkstra_time/1000.0);

    fprintf(stdout, "dijkstra-query: [query: %.2lfms %.2lfGiB/s] ", 
                     dijkstra_time, trans_rate/((double) (1<<30)));
    fprintf(stdout, "[trace: %.2lfms] [%.2lfms] ", trace_time, pop_time());
    print_pop_memtrack();
    fprintf(stdout, " ");
    print_current_mem();
    fprintf(stdout, "\n");
    fflush(stdout);
#endif

    FREE(d);
    FREE(visit);
    FREE(p);
    return min_cost;
}

/******************************************************* Program entry point. */

int main(int argc, char **argv)
{
#ifdef TRACK_RESOURCES
    push_time();
	push_memtrack();

    fprintf(stdout, "invoked as:");
    for(index_t f = 0; f < argc; f++) 
        fprintf(stdout, " %s", argv[f]);
    fprintf(stdout, "\n");
#endif

    if(argc > 1 && !strcmp(argv[1], "-h")) {
        fprintf(stdout, "Usage: %s -in <in-file> -src <source> -dst <destination>\n\n", argv[0]);
        return 0;
    }

    index_t has_input = 0;
    index_t has_source = 0;
    index_t has_destination = 0;
    char *filename = NULL;;
    index_t src = -1;
    index_t dst = -1;

    for(index_t f = 1; f < argc; f++) {
        if(argv[f][0] == '-') {
            if(!strcmp(argv[f], "-in")) {
                has_input = 1;
                filename = argv[++f];
            }
            if(!strcmp(argv[f], "-src")) {
                has_source = 1;
                src = atol(argv[++f]) - 1;
            }
            if(!strcmp(argv[f], "-dst")) {
                has_destination = 1;
                dst = atol(argv[++f]);
            }
        }
    }

    FILE *in = NULL;
    if(!has_input) {
        //read graph from standard input
#ifdef TRACK_RESOURCES
        fprintf(stdout, "Input file not specified, redirecting to standard input stream\n");
#endif
        in = stdin;
    } else {
        in = fopen(filename, "r");
        if(in == NULL)
            ERROR("unable to open file '%s'", filename);
    }

    if(!has_source) {
#ifdef TRACK_RESOURCES
        fprintf(stdout, "Source vertex not specified;");
        fprintf(stdout, " default source vertex is '0'\n");
        fflush(stdout);
#endif
        src = 0;
    }

    if(!has_destination) {
#ifdef TRACK_RESOURCES
        fprintf(stdout, "Destination vertex not specified;");
        fprintf(stdout, " reporting shortet path to all vertices.\n");
        fflush(stdout);
#endif
    } else {
        if(dst != -1) dst--;
    }

    //read input graph
    graph_t *g = graph_load(in);
    //build root query
    dijkstra_t *root = root_build(g, src, dst);
    //release graph memory
    graph_free(g);
    //execute the algorithm
    dijkstra_query(root);
    //release query memory
    dijkstra_free(root);

#ifdef TRACK_RESOURCES
    double time = pop_time();
    fprintf(stdout, "grand total [%.2lfms] ", time);
    print_pop_memtrack();
    fprintf(stdout, "\n");
    fprintf(stdout, "host: %s\n", sysdep_hostname());
    fprintf(stdout, "build: %s, %s\n",
                    "single thread",
                    "binary heap"
            );
    fprintf(stdout,
            "compiler: gcc %d.%d.%d\n",
            (index_t)__GNUC__,
            (index_t)__GNUC_MINOR__,
            (index_t)__GNUC_PATCHLEVEL__);
    fflush(stdout);
    assert(malloc_balance == 0);
    assert(memtrack_stack_top < 0);
#endif
    return 0;
}
