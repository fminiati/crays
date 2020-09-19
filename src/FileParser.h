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
#include <cstdlib>
#include <cstdlib>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

using std::cerr;
using std::endl;
using std::istringstream;
using std::map;
using std::string;

#ifndef _FILE_PARSER_H_
#define _FILE_PARSER_H_

class FileParser
{
public:
    // null constructor
    FileParser() : m_is_defined(false) {}

    // destructor
    ~FileParser() {}

    // full constructor
    FileParser(const string a_file_name)
        : m_is_defined(true)
    {
        // clear previous entries
        m_entries.clear();
        //
        std::ifstream file(a_file_name.c_str());
        if (!file)
        {
            std::cerr << "FileParser::define: Error while opening file!"
                      << '\n';
            exit(0);
        }

        string s;
        string key;
        while (getline(file, s))
        {
            // first non-blank char
            string::size_type pos = s.find_first_not_of(" ");

            // ignore comments
            if (pos != string::npos && s.at(pos) != '#')
            {
                pos = s.find("=");
                if (pos != string::npos)
                {
                    istringstream ikey(s.substr(0, pos));
                    string entry(s.substr(pos + 1, s.size() - 1));

                    // remove unwanted blanks
                    // string key;
                    ikey >> key;

                    m_entries[key] = entry;
                }
                else
                {
                    // add to previous key
                    m_entries[key].append(" " + s);
                }
            }
            else
            {
                // reset key
                key = " ";
            }
        }
        // close file
        file.close();
    }

    template <class T>
    void get_item(T &a_val, const string a_name)
    {
        assert(m_is_defined);

        // search object a_name
        auto it = m_entries.find(a_name);
        if (it != m_entries.end())
        {
            try
            {
                istringstream line(it->second);
                line >> a_val;
            }
            catch (std::ios_base::failure &)
            {
                std::cout << "FileParser::get_item: error: entry was : "
                          << it->second << endl;
            }
        }
        else
        {
            cerr << "FileParser::get_item: item " << a_name << " not found !" << endl;
            exit(0);
        }
    }

    // get a_n items into a_c container
    template <class T>
    void get_items(T &a_c, const int a_n, const string a_name)
    {
        assert(m_is_defined);

        // search object a_name
        auto it = m_entries.find(a_name.c_str());

        if (it != m_entries.end())
        {
            try
            {
                a_c.resize(a_n);

                istringstream line(it->second);
                for (auto &c : a_c)
                    line >> c; //*it;
            }
            catch (std::ios_base::failure &)
            {
                std::cout << "FileParser::get_item: error: entry was : "
                          << it->second << endl;
            }
        }
        else
        {
            cerr << "FileParser::get_item: item " << a_name << " not found !" << endl;
            exit(0);
        }

        if (a_n <= 0)
        {
            cerr << "FileParser::get_item: " << a_name
                 << " -- Warning: container size " << a_c.size() << endl;
        }
    }

protected:
    bool m_is_defined;
    map<string, string> m_entries;
};

#endif
