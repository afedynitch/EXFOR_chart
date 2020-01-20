from __future__ import print_function
import os
import sys
from os.path import join
import numpy as np
import periodictable
import six
from x4i3 import exfor_manager, exfor_entry
try:
    import cPickle as pickle
except:
    import pickle

data_dir = os.path.abspath('data')  # modified from above 

elem_names = {}
elem_Z = {}
for el in periodictable.elements:
    if el.number > 40:
        break
    elem_names[el.number] = el.symbol.upper()
    elem_Z[el.symbol.upper()] = el.number


def get_exfor_symb(nco_id):
    N, Z = get_N_Z(nco_id)
    return '{0}-{1}'.format(elem_names[Z], Z + N)


def get_N_Z(nco_id):
    Z = nco_id % 100
    A = (nco_id - Z)/100
    return A-Z, Z


def convert_to_full_tuple(nco_list, lim=30):
    outl = []
    for nco_id in nco_list:
        N, Z = get_N_Z(nco_id)
        if N > lim or Z > lim:
            continue
        name = '{0}-{1}'.format(elem_names[Z], Z+N)
        outl.append((N+Z, Z, nco_id, name))
    return outl


class NeucosDatabaseManager(object):

    def __init__(self, *args, **kwargs):
        self.endf = EndlDB('ENDF-B-VII.1', **kwargs)
        self.jendl = EndlDB('JENDL-PD-2004', **kwargs)
        # self.tendl = EndlDB('TENDL-2015', **kwargs)
        self.neucos = CombinedNeucosmaDB('TALYS_18_s1_gdr_patch',
                                         **kwargs)
        self.talys = NeucosmaDB('TALYS_18_default',**kwargs)
        self.psb = NeucosmaDB('PSB', **kwargs)
        self.peanut = NeucosmaDB('peanut', **kwargs)
        self.peanut_ias = NeucosmaDB('peanut_IAS', **kwargs)
        self.initialize()

    def initialize(self):
        self.endf_abs_elems = list(zip(*list(self.endf.db)))[1]
        self.jendl_abs_elems = list(zip(*list(self.jendl.db)))[1]
        # self.tendl_abs_elems = list(zip(*list(self.tendl.db)))[1]
        self.neucos_abs_elems = [el[1] for el in self.neucos.db
                                 if el[0] == 0]
        self.talys_abs_elems = [el[1] for el in self.talys.db
                                 if el[0] == 0]
        self.psb_abs_elems = [el[1] for el in self.psb.db
                              if el[0] == 0]
        self.peanut_abs_elems = [el[1] for el in self.peanut.db
                              if el[0] == 0]
        self.peanut_ias_abs_elems = [el[1] for el in self.peanut_ias.db
                              if el[0] == 0]                   


    def intersection(self, db1_elems, db2_elems):
        return list(set(db1_elems).intersection(
            db2_elems))


class ExforDatabaseManager(object):

    def __init__(self, debug=False, **kwargs):
        print('Loading database x4i3 Database manager.')
        self.x4man = exfor_manager.X4DBManagerPlainFS()
        self.d = debug
        self.reaction = 'G,*'
        self.observable = 'SIG'
        self.filters = dict(EVAL=False, PAR=False,
                            EP=False, UKNTAR=False,
                            COMBI=False)
        self.select_x4keys = {}
        self.select_ncokeys = {}

    def set_reaction(self, cstype):
        self.reaction = 'G,' + cstype

    def set_filter(self, filter, value):
        self.filters[filter] = value

    def get_entries(self, nco_id):
        target = get_exfor_symb(nco_id) if nco_id > 0 else nco_id
        query = dict(
            target=target,
            reaction=self.reaction,
            quantity=self.observable)
        if nco_id < 0:
            query.pop('target')

        subents = self.x4man.retrieve(**query)
        
        self.select_x4keys = {}
        self.select_ncokeys = {}
        for key, x4ent in six.iteritems(subents):
            ds = x4ent.getDataSets()
            for key, dstr in six.iteritems(ds):
                kred = key[:-1]
                if kred in self.select_x4keys:
                    continue
                e = EXFOREntry(kred, dstr, debug=self.d)
        #         print(e)
                if not e.is_useful(**self.filters):
                    continue
                self.select_ncokeys[e.get_db_input()] = e
                self.select_x4keys[kred] = e
        print('Found {0} entries'.format(len(self.select_x4keys)))
        return self.select_x4keys

    def plot(self, index, on_err=False, **kwargs):
        import matplotlib.pyplot as plt
        if None in index:
            return
        if index not in self.select_ncokeys:
            return
        self.select_ncokeys[index].plot(**kwargs)

    def get_nco_list_of_selection(self):
        return list(zip(*sorted(self.select_ncokeys)))[1]


class NuclearDB(object):

    # def plot_list(self, p_list, on_err=False, **kwargs):
    #     import matplotlib.pyplot as plt
    #     if None in p_list:
    #         return
    #     d = None
    #     for mid, did in p_list:
    #         egr, di = self.get_data(mid, did)
    #         if d == None:
    #             d = np.zeros_like(egr)
    #         if np.sum(di) > 0:
    #             d += di
    #         else:
    #             print('something wrong with data')
    #             return

    #     plt.plot(egr, d, **kwargs)

    def get_data(self, mid, did, on_err=False):
        if (mid, did) in self.db:
            return self.db[(mid, did)]
        else:
            if on_err:
                raise Exception(
                    'Mother/daughter combination {0}/{1} not found'.format(
                        mid, did))
            else:
                return np.zeros(5), np.zeros(5)

    def plot(self, index, on_err=False, ret_data=False, **kwargs):
        import matplotlib.pyplot as plt
        mid, did = index
        egr, d = self.get_data(mid, did, on_err)
        if sum(d) <= 0:
            return
        if 'axes' in kwargs:
            kwargs['axes'].plot(egr, d, **kwargs)
        else:
            plt.plot(egr, d, **kwargs)


class EndlDB(NuclearDB):

    def __init__(self, db_prefix, debug=False, **kwargs):
        self.db = {}
        self.d = debug
        self.prefix = db_prefix
        self.load_and_parse()

    def load_and_parse(self):
        try:
            self.db = pickle.load(open(
                join(data_dir, self.prefix + '-db.ppd'), 'rb'), encoding='latin1')
        except:
            self.db = pickle.load(open(
                join(data_dir, self.prefix + '-db.ppd'), 'rb'))

    def plot(self, index, on_err=False, **kwargs):
        # Only plots first element, no sums
        import matplotlib.pyplot as plt
        if None in index:
            return
        egr, d = self.get_data(*index)
        if sum(d) <= 0.:
            return
        if 'axes' in kwargs:
            kwargs['axes'].plot(egr, d, **kwargs)
        else:
            plt.plot(egr, d, **kwargs)

class NeucosmaDB(NuclearDB):

    def __init__(self, db_prefix, debug=False, **kwargs):
        self.db = {}
        self.d = debug
        self.db_prefix = db_prefix
        self.load_or_parse()

    def load_or_parse(self):
        prefix = self.db_prefix

        try:
            self.db = pickle.load(open(join(data_dir, prefix + '_db.ppd'), 'rb'))
            self.zeros = np.zeros_like(list(self.db.items())[0][0])
            return
        except UnicodeDecodeError:
            self.db = pickle.load(open(join(data_dir, prefix + '_db.ppd'), 'rb'),
                encoding='latin1')
            self.zeros = np.zeros_like(list(self.db.items())[0][0])
            return
        except IOError:
            print('{0}::load_or_parse(): no precompiled DB found. Parsing ASCII...'
                .format(self.__class__.__name__))

        # If no pre-compiled database found, parse the ASCII files

        f_nonel = join(data_dir, prefix + '_nonel.dat')
        f_incl = join(data_dir, prefix + '_incl_i_j.dat')
        egrid_fname = join(data_dir, prefix.split('_')[0] +
                           '_e-grid.grid')

        for fn in [egrid_fname, f_incl, f_nonel]:
            if not os.path.isfile(fn):
                raise Exception('File {0} not found.'.format(fn))

        self.e_grid = np.loadtxt(egrid_fname)
        self.zeros = np.zeros_like(self.e_grid)
        nonel = np.loadtxt(f_nonel)
        self.db = {}
        for lidx in range(nonel.shape[0]):
            self.db[(0, int(nonel[lidx, 0]))] = (self.e_grid, nonel[lidx, 1:])

        incl = np.loadtxt(f_incl)
        for lidx in range(incl.shape[0]):
            mid = int(incl[lidx, 0])
            did = int(incl[lidx, 1])
            self.db[(mid, did)] = (self.e_grid, incl[lidx, 2:])

        # Dump databse to file
        pickle.dump(self.db,
                    open(join(data_dir, prefix + '_db.ppd'), 'wb'),
                    protocol=-1)


class CombinedNeucosmaDB(NeucosmaDB):

    def __init__(self, db_prefix):
        self.talys_db = NeucosmaDB(db_prefix)
        self.crp2_db = NeucosmaDB('CRP2')
        self.db_prefix = db_prefix + '+' 'CRP2'
        self.db = self.talys_db.db
        self.db.update(self.crp2_db.db)
        del self.talys_db, self.crp2_db


class ReactionCombination(object):

    def __init__(self, reac_list):
        self.reactions = [Reaction(r) for r in reac_list]
        self.targ = self.reactions[0].targ
        self.comb_str = str(reac_list.__repr__())

    def __repr__(self):
        rstr = 'Reaction list:' + self.comb_str + '\n'
        for i, r in enumerate(self.reactions):
            rstr += str(i) + ":, " + '\n'
            rstr += r.__repr__(o=1) + '\n'
        return rstr

    def is_useful(self, **kwargs):
        for r in self.reactions:
            if not r.is_useful(**kwargs):
                return False
        return True


class Reaction(object):

    def __init__(self, x4reac):
        self.quantities = x4reac.quantity
        self.proj = x4reac.proj
        self.targ = self.prod2nco_id(x4reac.targ)
        self.prod_list = [self.prod2nco_id(p) for p in x4reac.products]
        self.proc_type = x4reac.processType
        self.residual = self.prod2nco_id(x4reac.residual)
        self.quantities = x4reac.quantity
        self.r_str = str(x4reac.__repr__())

    def prod2nco_id(self, p):
        try:
            return p.getA() * 100 + p.getZ()
        except:
            return p

    def is_useful(self, **kwargs):
        for key, value in six.iteritems(kwargs):
            if not value and key in self.quantities:
                print('Rejecting due to filter: ', key)
                return False
            if not value and key == 'UKNTAR' and \
                    abs(self.targ) > 7000:
                print('Rejecting due to UKNTAR: ', self.targ)
                return False
        return True

    def __repr__(self, o=0):
        offs = o*'   '
        rstr = ''
        rstr += offs + 'Reaction str:' + str(self.r_str) + '\n'
        rstr += offs + 'projectile  :' + str(self.proj) + '\n'
        rstr += offs + 'target      :' + str(self.targ) + '\n'
        rstr += offs + 'process type:' + str(self.proc_type) + '\n'
        rstr += offs + 'products    :' + str(self.prod_list) + '\n'
        for p in self.prod_list:
            rstr += offs + '         -> ' + str(p) + '\n'
        rstr += offs + 'residual    :' + str(self.residual) + '\n'
        rstr += offs + 'quantities  :' + str(self.quantities) + '\n'
        return rstr


class EXFOREntry(object):

    def __init__(self, key, dataset, debug=False):
        self.entry_key, self.subentry_key = key
        r = dataset.reaction.__str__()
        if hasattr(dataset.reaction[0], 'reaction_list'):
            self.reaction = ReactionCombination(
                dataset.reaction[0].reaction_list)
            self.is_sum = True
        else:
            self.reaction = Reaction(dataset.reaction[0])
            self.is_sum = False

        self.d = debug
        self.parse_success = self.parse(dataset)

        self.blacklist = ['M0372005','M0372004']

    def is_useful(self, **kwargs):
        if self.subentry_key in self.blacklist:
            print('Dataset blacklisted', self.subentry_key)
            return False
        if not self.parse_success:
            if self.d > 1:
                print('Could not parse data for', self.subentry_key)
            return False
        if 'COMBI' in kwargs and not kwargs['COMBI'] and self.is_sum:
            print('Rejecting reaction combination')
            return False

        return self.reaction.is_useful(**kwargs)

    def parse(self, dataset):

        if self.d > 2:
            print('Parsing ', self.entry_key, self.subentry_key)
        cols = dataset.numcols()
        units = [u.upper() for u in dataset.units]
        labels = [l.upper() for l in dataset.labels]

        self.entr_label = ('{0} et al.,{1},{2}'.format(
            dataset.author[0].split('.')[-1],
            dataset.year, self.subentry_key[-1]))

        self.publabel = ('{0} et al.,{1}'.format(
            dataset.author[0].split('.')[-1],
            dataset.year))

        if self.d > 1:
            print('\t cols:', dataset.numcols())
            print('\t units:', dataset.units)
            print('\t reac:', dataset.reaction)
            print('\t labels:', dataset.labels)

        en_col = 0  # column for energy
        data_col = 0  # column for data
        err_col = 0  # column for error (not implmemented)

        conv_e = 1.  # conversion factor to MeV
        conv_cs = 1.  # conversion factor to mb

        if 'EN' in labels:
            en_col = labels.index('EN')
        elif 'EN-MAX' in labels:
            if self.d > 1:
                print('Using EN-MAX')
            en_col = labels.index('EN-MAX')
        else:
            if self.d > 1:
                print('Energy not found in ' +
                  ','.join(labels))
                print(dataset)
            return False

        if 'DATA' in labels:
            data_col = labels.index('DATA')
        elif 'DATA-MAX' in labels:
            print('Using DATA-MAX')
            data_col = labels.index('DATA-MAX')
        else:
            print('DATA not found in ' +
                  ','.join(labels))

            if self.d:
                print(self)
            return False

        if units[en_col] == 'GEV':
            conv_e = 1e3
        elif units[en_col] == 'KEV':
            conv_e = 1e-3
        elif units[en_col] != 'MEV':
            if self.d > 0:
                print(dataset)
            raise Exception(
                'Unknown energy unit ' + units[en_col])
        if units[data_col] == 'BARNS':
            conv_cs = 1e3
        elif units[data_col] == 'MICRO-B':
            conv_cs = 1e-3
        elif units[data_col] != 'MB':
            if self.d > 1:
                print('Unknown data unit ' + units[data_col])
                print(dataset)
            return False

        dtable = np.array(list(zip(*dataset.data)), dtype=np.float)

        # ma_data = np.ma.mask_invalid(dtable[data_col])
        # msk = ma_data.mask
        self.e_grid = conv_e*dtable[en_col]
        self.data = conv_cs*dtable[data_col]

        if self.d > 3:
            print(dataset)

        return True

    def get_db_input(self):
        if isinstance(self.reaction, ReactionCombination):
            return [(r.targ, r.residual)
                    for r in self.reaction.reactions]
        else:
            if self.reaction.residual != None:
                return (self.reaction.targ, self.reaction.residual)
            elif (self.reaction.prod_list[0] == 'X' and
                    self.reaction.prod_list[1] > 0):
                print('X+xx found', self.reaction.prod_list[1])
                print(self.reaction)
                return (self.reaction.targ, self.reaction.prod_list[1])
            elif self.reaction.prod_list[0] == 'ABS':
                return (0, self.reaction.targ)
            else:
                print('No residual defined')
                print(self.reaction)
                return None

    def plot(self, label='author', **kwargs):
        import matplotlib.pyplot as plt
        if label == 'short':
            kwargs['label'] = self.publabel
        elif isinstance(self.reaction, ReactionCombination):
            kwargs['label'] = self.reaction.comb_str
        else:
            kwargs['label'] = self.reaction.r_str
        if 'axes' in kwargs:
            kwargs.pop('axes').errorbar(self.e_grid, self.data, 
                ms=3.5, **kwargs)
        else:
            plt.errorbar(self.e_grid, self.data, ms=3.5, **kwargs)
