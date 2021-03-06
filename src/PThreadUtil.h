//
// Copyright (C) 2020 Francesco Miniati <francesco.miniati@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
#ifndef PTHREADUTIL_H
#define PTHREADUTIL_H

#include <vector>
#include <pthread.h>
#include <functional>
#include <cassert>

namespace fm::cr_transport
{
    template <typename Op, typename Obj, typename V>
    struct ThreadUt
    {
        size_t _num_threads;
        size_t _loop_size;

        ThreadUt() = default;
        ThreadUt(ThreadUt const &) = default;

        // template thread arg struct
        struct ThreadArg
        {
            ThreadArg() = default;
            ThreadArg(ThreadArg const &) = default;

            size_t beg, strd, size;
            Op *op;
            Obj *obj;
            V *var;
        };

        /// aliases
        using rw_obj_t = std::reference_wrapper<Obj>;
        using rw_var_t = std::reference_wrapper<V>;

        size_t launch_threads(const Op &a_op, Obj &a_obj, std::vector<V> &a_v)
        {
            std::vector<rw_var_t> v(a_v.begin(), a_v.end());
            std::vector<rw_obj_t> obj(_num_threads, std::ref(a_obj));
            return launch_threads(a_op, obj, v);
        }

        ///
        size_t launch_threads(const Op &a_op, std::vector<Obj> &a_obj, V &a_v)
        {
            std::vector<rw_var_t> v(_num_threads, std::ref(a_v));
            std::vector<rw_obj_t> obj(a_obj.begin(), a_obj.end());
            return launch_threads(a_op, obj, v);
        }

        ///
        size_t launch_threads(const Op &a_op, Obj &a_obj, V &a_v)
        {
            std::vector<rw_var_t> v(_num_threads, std::ref(a_v));
            std::vector<rw_obj_t> obj(_num_threads, std::ref(a_obj));
            return launch_threads(a_op, obj, v);
        }

        ///
        size_t launch_threads(const Op &a_op, std::vector<rw_obj_t> &a_obj, std::vector<rw_var_t> &a_v)
        {
            pthread_attr_t th_attr;
            pthread_attr_init(&th_attr);
            pthread_attr_setdetachstate(&th_attr, PTHREAD_CREATE_JOINABLE);

            // create threads
            std::vector<pthread_t> th(_num_threads);
            std::vector<ThreadArg> ta(_num_threads);

            auto t_lambda = [](void *args) -> void * {
                ThreadArg *ta = static_cast<ThreadArg *>(args);
                ta->op->operator()(*ta->obj, *ta->var, ta->beg, ta->strd, ta->size);
                return nullptr;
            };

            size_t err = 0;
            for (size_t t = 0; t < _num_threads; ++t)
            {
                ta[t].beg = t;
                ta[t].strd = _num_threads;
                ta[t].size = _loop_size;
                ta[t].op = new Op(a_op);
                ta[t].obj = &a_obj[t].get();
                ta[t].var = &a_v[t].get();
                err += pthread_create(&th[t], &th_attr, t_lambda, (void *)&ta[t]);
            }

            void *status;
            for (auto &t : th)
            {
                err += pthread_join(t, &status);
            }
            for (auto &a : ta)
            {
                delete a.op;
                a.op = nullptr;
            }
            assert(err == 0);
            return err;
        }
    };

}; // namespace fm::cr_transport

#endif
