import { Injectable } from '@angular/core';
import { Firestore } from '@angular/fire/firestore';
import { User } from '../models/user';

@Injectable({
  providedIn: 'root'
})
export class UserService {

  constructor(private firestore:Firestore) { }

  createUser(user: User) {
    // create user
  }
}
