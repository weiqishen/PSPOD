/**
 * @file ndarray.h
 * @brief A simple N-dimensional array in Column Major
 * @author Weiqi Shen weiqishen1994@ufl.edu
 * @version 0.1
 * @date 2019-01-26
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#pragma once
#include <algorithm>
#include <initializer_list>
#include <iostream>
#include "err.h"

/**
 * @brief A simple N-dimensional array in Column Major
 * @class ndarray
 * @tparam T data type
 */
template <typename T>
class ndarray
{
public:
  // #### constructors ####

  /// Default constructor to construct an empty ndarray
  ndarray(void);

  /**
   * @brief N-Dimension constructor
   * 
   * @param list A list to specify dimension of the ndarray
   */
  ndarray(std::initializer_list<int> list);

  /**
   * @brief 1-D constructor
   * 
   * @param nele number of element
   */
  ndarray(int nele);

  /// Copy constructor
  ndarray(const ndarray<T> &in_array);

  /// destructor
  ~ndarray(void);

  // #### methods ####

  /**
   * @brief Overload stream output operator
   * 
   * @param out ostream object
   * @param s array to output
   * @return ostream& reference to the ostream object
   */
  template<class U>
  friend std::ostream& operator<<(std::ostream& out, ndarray<U>& s);

  /**
   * @brief setup a new ndarray object
   * 
   * @param list A list to specify dimension of the ndarray
   */
  void setup(std::initializer_list<int> list);

  /**
   * @brief setup a new 1-D ndarray object
   * 
   * @param list A list to specify dimension of the ndarray
   */
  void setup(int nele);

  ///Assignment
  ndarray<T> &operator=(const ndarray<T> &in_array);

  ///fill with value
  ndarray<T> &operator=(const T &other);

  /// Access/set ndarray element
  T &operator()(std::initializer_list<int> list);

  /// 1-D access/set
  T &operator()(int idx);

  /// Return pointer of ndarray element
  T *get_ptr(std::initializer_list<int> list);

  /// Return pointer of ndarray element in 1-D
  T *get_ptr(int idx = 0);

  /// Get number of elements along one axis
  int get_dim(int in_dim);

  /// Get number of dimension of the ndarray
  int get_n_dim(void);

  ///get the length of the ndarray
  int get_len(void);

  /// Method to get maximum value of ndarray
  T get_max(void);

  /// Method to get minimum value of ndarray
  T get_min(void);

  /// Reshape the array
  void reshape(std::initializer_list<int> list);

protected:
  int *shape;
  T *data;
  int n_dim;
  int len;

private:
  /// helper method to calculate length of ndarray
  void calc_len(void);
};

// definitions

using namespace std;

// #### constructors ####

// default constructor

template <typename T>
ndarray<T>::ndarray()
{
  len = 0;
  n_dim = 0;
  shape = NULL;
  data = NULL;
}

// constructor 1

template <typename T>
ndarray<T>::ndarray(initializer_list<int> list)
{
  //store dimension array
  n_dim = list.size();
  shape = new int[n_dim];

  int i = 0;
  for (auto l : list)
    shape[i++] = l;

  calc_len();
  data = new T[len];
}

// constructor 2

template <typename T>
ndarray<T>::ndarray(int nele)
{
  //store dimension array
  n_dim = 1;
  shape = new int[n_dim];
  shape[0] = nele;
  calc_len();
  data = new T[len];
}

// copy constructor

template <typename T>
ndarray<T>::ndarray(const ndarray<T> &in_array)
{
  n_dim = in_array.n_dim;
  shape = new int[n_dim];
  copy(in_array.shape, in_array.shape + n_dim, this->shape);

  calc_len();
  data = new T[len];
  copy(in_array.data, in_array.data + len, this->data);
}

// assignment

template <typename T>
ndarray<T> &ndarray<T>::operator=(const ndarray<T> &in_array)
{

  if (this == &in_array)
  {
    return (*this);
  }
  else
  {
    delete[] data;
    delete[] shape;

    n_dim = in_array.n_dim;
    shape = new int[n_dim];
    copy(in_array.shape, in_array.shape + n_dim, this->shape);

    calc_len();
    data = new T[len];
    copy(in_array.data, in_array.data + len, this->data);
    return (*this);
  }
}

template <typename T>
ndarray<T> &ndarray<T>::operator=(const T &other)
{
  fill_n(this->data, len, other);
}
// destructor

template <typename T>
ndarray<T>::~ndarray()
{
  delete[] data;
  delete[] shape;
}

// #### methods ####

// setup

template <typename T>
void ndarray<T>::setup(initializer_list<int> list)
{
  //delete previous data
  delete[] shape;
  delete[] data;

  //store dimension array
  n_dim = list.size();
  shape = new int[n_dim];

  int i = 0;
  for (auto l : list)
    shape[i++] = l;

  calc_len();
  data = new T[len];
}

//1-D setup
template <typename T>
void ndarray<T>::setup(int nele)
{
  //delete previous data
  delete[] shape;
  delete[] data;

  //store dimension array
  n_dim = 1;
  shape = new int[n_dim];
  shape[0]=nele;
  calc_len();
  data = new T[len];
}

template <typename T>
T &ndarray<T>::operator()(int idx)
{
  return data[idx];
}

template <typename T>
T &ndarray<T>::operator()(initializer_list<int> list)
{
  int idx = 0, acc = 1;
  int i = 0;

#ifdef _DEBUG
  if (list.size() != n_dim)
    Fatal_Error("Invalid dimension");
#endif

        for (auto l : list)
    {
#ifdef _DEBUG
      if (l >= shape[i])
        Fatal_Error("Out of bound");
#endif
      idx += acc * l;
      acc *= shape[i++];
    }

  return data[idx];
}

// return pointer

template <typename T>
T *ndarray<T>::get_ptr(int idx)
{
  return data + idx;
}

template <typename T>
T *ndarray<T>::get_ptr(initializer_list<int> list)
{
  int idx = 0, acc = 1;
  int i = 0;

#ifdef _DEBUG
  if (list.size() != n_dim)
    Fatal_Error("Invalid dimension");
#endif

        for (auto l : list)
    {
#ifdef _DEBUG
      if (l >= shape[i])
        Fatal_Error("Out of bound");
#endif
      idx += acc * l;
      acc *= shape[i++];
    }

  return data + idx;
}

// obtain dimension

template <typename T>
int ndarray<T>::get_dim(int in_dim)
{
  if (in_dim < n_dim)
    return shape[in_dim];
  else
    Fatal_Error("Dimension not supported");
}

template <typename T>
int ndarray<T>::get_n_dim()
{
  return n_dim;
}

template <typename T>
void ndarray<T>::calc_len(void)
{
  len = 1;
  for (int i = 0; i < n_dim; i++)
    len *= shape[i];
}

template <typename T>
int ndarray<T>::get_len(void)
{
  return len;
}
// method to calculate maximum value of ndarray
// Template specialization
template <typename T>
T ndarray<T>::get_max(void)
{
  return *max_element(data, data + len);
}

// method to calculate minimum value of ndarray
// Template specialization
template <typename T>
T ndarray<T>::get_min(void)
{
  return *min_element(data, data + len);
}

template <typename T>
void ndarray<T>::reshape(initializer_list<int> list)
{

  int i = 0;

  if (n_dim != list.size())
  {
    n_dim = list.size();
    //delete shape
    delete[] shape;
    shape=new int[n_dim];
  }
  
#ifdef _DEBUG
  int acc = 1;
#endif

  for (auto l : list)
  {
#ifdef _DEBUG
    acc *= l;
#endif
    shape[i++] = l;
  }

#ifdef _DEBUG
  if (acc != len)
    Fatal_Error("Total number of element doesn't agree");
#endif
}

//-------------friend---------------------
//output
template <typename U>
ostream &operator<<(ostream &out, ndarray<U> &s)
{
  if (s.n_dim == 1) //1-D output
  {
    for (int i = 0; i < s.shape[0]; i++)
      out << setw(10) << s(i);
    out << endl;
  }
  else if (s.n_dim == 2) //2d output
  {
    for (int i = 0; i < s.shape[0]; i++)
    {
      for (int j = 0; j < s.shape[1]; j++)
      {
        out << setw(10) << s({i, j});
      }
      out << endl;
    }
  }
  else if (s.n_dim == 3) //3d output
  {
    for (int k = 0; k < s.shape[2]; k++)
    {
      out << "slice: (:,:," << k << ")" << endl;
      for (int i = 0; i < s.shape[0]; i++)
      {
        for (int j = 0; j < s.shape[1]; j++)
        {
          out << setw(10) << s({i, j, k});
        }
        out << endl;
      }
    }
  }
  else
  {
    Fatal_Error("Output not supported");
  }

  return out;
}